#!/bin/bash
#SBATCH --job-name=ANI_dRep
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --mem=128GB
#SBATCH --output=/path/to/crc-metagenome-mag-pipeline/C_genomes/slurm_out/slurm_ANI_dRep.out
#SBATCH --error=/path/to/crc-metagenome-mag-pipeline/C_genomes/slurm_out/slurm_ANI_dRep.err

set -euo pipefail

###############################################################################
# 方案B：C1A + C1B 全部合并，用 dRep 按 99% ANI 去冗余
# 标签规则
#   - 只要同一个 secondary_cluster 中同时出现 C1A 和 C1B -> cluster_type = MIXED
#   - 否则只有 C1A -> cluster_type = C1A
#   - 否则只有 C1B -> cluster_type = C1B
#
# 最终输出三类 winners 代表基因组到三个文件夹：
#   final_reps/C1A
#   final_reps/C1B
#   final_reps/MIXED
###############################################################################

# ===== Conda 环境=====
CONDA_BIN="/path/to/conda/condabin/conda"
eval "$("$CONDA_BIN" shell.bash hook)"
set +u
conda activate CRC
set -u

# ====== 路径与参数 ======
BASE_DIR="/path/to/crc-metagenome-mag-pipeline/C_genomes"
C1A_DIR="${BASE_DIR}/C1A_genomes"
C1B_DIR="${BASE_DIR}/C1B_genomes"

SA="0.99"           # secondary ANI threshold
PA="0.90"           # primary Mash threshold
NC="0.10"           # minimum alignment coverage
THREADS="128"       # 并行线程

# 强制不跑 CheckM
IGNORE_GENOME_QUALITY="true"

###############################################################################
# 0) 基本检查
###############################################################################
command -v dRep >/dev/null 2>&1 || { echo "[ERROR] dRep not found in PATH"; exit 1; }
command -v python3 >/dev/null 2>&1 || { echo "[ERROR] python3 not found in PATH"; exit 1; }

[[ -d "${C1A_DIR}" ]] || { echo "[ERROR] C1A_DIR not found: ${C1A_DIR}"; exit 1; }
[[ -d "${C1B_DIR}" ]] || { echo "[ERROR] C1B_DIR not found: ${C1B_DIR}"; exit 1; }

TS="$(date +%Y%m%d_%H%M%S)"
WORKDIR="${BASE_DIR}/drep_all_sa${SA}_${TS}"
LINKDIR="${WORKDIR}/ALL_genomes_links"
DREP_OUT="${WORKDIR}/dRep_out"
REPORT_DIR="${WORKDIR}/reports"
FINAL_REP_DIR="${WORKDIR}/final_reps"

mkdir -p "${WORKDIR}" "${LINKDIR}" "${REPORT_DIR}" "${FINAL_REP_DIR}/C1A" "${FINAL_REP_DIR}/C1B" "${FINAL_REP_DIR}/MIXED"

echo "[INFO] WORKDIR      = ${WORKDIR}"
echo "[INFO] LINKDIR      = ${LINKDIR}"
echo "[INFO] DREP_OUT     = ${DREP_OUT}"
echo "[INFO] REPORT_DIR   = ${REPORT_DIR}"
echo "[INFO] FINAL_REPS   = ${FINAL_REP_DIR}"

###############################################################################
# 1) 收集基因组并建立“带标签前缀”的软链接，避免同名冲突
#    同时生成 labels.tsv：genome<tab>orig_cluster<tab>orig_path
###############################################################################
LABELS_TSV="${WORKDIR}/labels.tsv"
echo -e "genome\torig_cluster\torig_path" > "${LABELS_TSV}"

shopt -s nullglob

C1A_FILES=("${C1A_DIR}"/*.fa "${C1A_DIR}"/*.fna "${C1A_DIR}"/*.fasta "${C1A_DIR}"/*.fas)
C1B_FILES=("${C1B_DIR}"/*.fa "${C1B_DIR}"/*.fna "${C1B_DIR}"/*.fasta "${C1B_DIR}"/*.fas)

if [[ ${#C1A_FILES[@]} -eq 0 ]]; then
  echo "[ERROR] No genome fasta found in ${C1A_DIR}"
  exit 1
fi
if [[ ${#C1B_FILES[@]} -eq 0 ]]; then
  echo "[ERROR] No genome fasta found in ${C1B_DIR}"
  exit 1
fi

echo "[INFO] Found C1A genomes: ${#C1A_FILES[@]}"
echo "[INFO] Found C1B genomes: ${#C1B_FILES[@]}"

link_one () {
  local orig="$1"
  local label="$2"   # C1A or C1B
  local base
  base="$(basename "$orig")"
  local new="${label}__${base}"
  ln -sf "$orig" "${LINKDIR}/${new}"
  echo -e "${new}\t${label}\t${orig}" >> "${LABELS_TSV}"
}

for f in "${C1A_FILES[@]}"; do
  link_one "$f" "C1A"
done

for f in "${C1B_FILES[@]}"; do
  link_one "$f" "C1B"
done

echo "[INFO] labels.tsv written: ${LABELS_TSV}"

###############################################################################
# 2) 运行 dRep 99% ANI 去冗余
###############################################################################
GENOME_LINKS=("${LINKDIR}"/*.fa "${LINKDIR}"/*.fna "${LINKDIR}"/*.fasta "${LINKDIR}"/*.fas)
if [[ ${#GENOME_LINKS[@]} -eq 0 ]]; then
  echo "[ERROR] No linked genomes found in ${LINKDIR}"
  exit 1
fi

echo "[INFO] Running dRep dereplicate (NO CheckM)..."

DREP_CMD=(dRep dereplicate "${DREP_OUT}"
  -g "${GENOME_LINKS[@]}"
  -sa "${SA}"
  -pa "${PA}"
  -nc "${NC}"
  -p "${THREADS}"
  --S_algorithm fastANI
  --ignoreGenomeQuality
)

echo "[INFO] dRep command:"
printf "  %q" "${DREP_CMD[@]}"
echo

"${DREP_CMD[@]}"

echo "[INFO] dRep finished."

###############################################################################
# 3) 解析 dRep 输出（Cdb/Wdb）+ labels.tsv，生成报告（不使用 majority）
###############################################################################
CDB="${DREP_OUT}/data_tables/Cdb.csv"
WDB="${DREP_OUT}/data_tables/Wdb.csv"

if [[ ! -s "${CDB}" ]]; then
  echo "[ERROR] Cdb.csv not found or empty: ${CDB}"
  exit 1
fi
if [[ ! -s "${WDB}" ]]; then
  echo "[ERROR] Wdb.csv not found or empty: ${WDB}"
  exit 1
fi

PY="${WORKDIR}/make_reports_and_split3.py"

cat > "${PY}" << 'PYCODE'
#!/usr/bin/env python3
import os
import sys
import pandas as pd

def die(msg: str, code: int = 1):
    print(f"[ERROR] {msg}", file=sys.stderr)
    sys.exit(code)

def classify_cluster(c1a: int, c1b: int) -> str:
    # Classification rule: mixed membership takes precedence over C1A or C1B alone.
    if c1a > 0 and c1b > 0:
        return "MIXED"
    if c1a > 0 and c1b == 0:
        return "C1A"
    if c1b > 0 and c1a == 0:
        return "C1B"
    return "UNKNOWN"

def main():
    cdb = os.environ.get("CDB")
    wdb = os.environ.get("WDB")
    labels_tsv = os.environ.get("LABELS_TSV")
    report_dir = os.environ.get("REPORT_DIR")
    drep_rep_dir = os.environ.get("DREP_REP_GENOMES_DIR")
    final_rep_dir = os.environ.get("FINAL_REP_DIR")

    for fp in [cdb, wdb, labels_tsv]:
        if not fp or not os.path.exists(fp):
            die(f"Missing file: {fp}")

    os.makedirs(report_dir, exist_ok=True)
    os.makedirs(final_rep_dir, exist_ok=True)

    # labels.tsv: genome orig_cluster orig_path
    lab = pd.read_csv(labels_tsv, sep="\t")
    if not {"genome", "orig_cluster"}.issubset(set(lab.columns)):
        die("labels.tsv must include columns: genome, orig_cluster")

    # Cdb.csv
    c = pd.read_csv(cdb)
    if not {"genome", "secondary_cluster"}.issubset(set(c.columns)):
        die("Cdb.csv must include columns: genome, secondary_cluster")

    # Wdb.csv
    w = pd.read_csv(wdb)
    if not {"genome", "cluster"}.issubset(set(w.columns)):
        die("Wdb.csv must include columns: genome, cluster")

    # Merge labels into Cdb
    df = c.merge(lab[["genome", "orig_cluster"]], on="genome", how="left")
    df["orig_cluster"] = df["orig_cluster"].fillna("UNKNOWN")

    # Cluster composition
    comp = (
        df.groupby("secondary_cluster")["orig_cluster"]
        .value_counts()
        .unstack(fill_value=0)
        .reset_index()
    )
    for col in ["C1A", "C1B"]:
        if col not in comp.columns:
            comp[col] = 0

    comp_size = df.groupby("secondary_cluster").size().reset_index(name="n_genomes")
    comp = comp.merge(comp_size, on="secondary_cluster", how="left")

    comp["cluster_type"] = comp.apply(lambda r: classify_cluster(int(r["C1A"]), int(r["C1B"])), axis=1)
    comp["mixed_C1A_C1B"] = comp["cluster_type"].eq("MIXED")

    # Winners mapping: Wdb.cluster == secondary_cluster
    winners = w.rename(columns={"cluster": "secondary_cluster", "genome": "winner_genome"})
    winners = winners[["secondary_cluster", "winner_genome"]].drop_duplicates()

    # Original C1A/C1B label, retained for provenance but not used for classification.
    winners = winners.merge(
        lab[["genome", "orig_cluster"]].rename(columns={"genome": "winner_genome", "orig_cluster": "winner_orig_label"}),
        on="winner_genome", how="left"
    )
    winners["winner_orig_label"] = winners["winner_orig_label"].fillna("UNKNOWN")

    # Cluster summary
    cluster_summary = comp.merge(winners, on="secondary_cluster", how="left")
    cluster_summary = cluster_summary[
        ["secondary_cluster", "n_genomes", "C1A", "C1B", "cluster_type", "mixed_C1A_C1B", "winner_genome", "winner_orig_label"]
    ].sort_values(["mixed_C1A_C1B", "n_genomes"], ascending=[False, False])

    cluster_summary_fp = os.path.join(report_dir, "cluster_summary.tsv")
    cluster_summary.to_csv(cluster_summary_fp, sep="\t", index=False)

    # Genome to cluster map
    df2 = df.merge(comp[["secondary_cluster", "C1A", "C1B", "cluster_type", "mixed_C1A_C1B"]], on="secondary_cluster", how="left")
    genome_cluster_map = df2[[
        "genome", "orig_cluster", "secondary_cluster",
        "cluster_type", "mixed_C1A_C1B", "C1A", "C1B"
    ]].sort_values(["mixed_C1A_C1B", "secondary_cluster"], ascending=[False, True])

    genome_cluster_map_fp = os.path.join(report_dir, "genome_cluster_map.tsv")
    genome_cluster_map.to_csv(genome_cluster_map_fp, sep="\t", index=False)

    # Mixed cluster genomes list
    mixed_cluster_genomes = genome_cluster_map[genome_cluster_map["cluster_type"].eq("MIXED")]
    mixed_cluster_genomes_fp = os.path.join(report_dir, "mixed_cluster_genomes.tsv")
    mixed_cluster_genomes.to_csv(mixed_cluster_genomes_fp, sep="\t", index=False)

    # Winners with cluster type
    winners_with_type = cluster_summary[[
        "secondary_cluster", "cluster_type", "winner_genome", "winner_orig_label", "n_genomes", "C1A", "C1B", "mixed_C1A_C1B"
    ]]
    winners_with_type_fp = os.path.join(report_dir, "winners_with_cluster_type.tsv")
    winners_with_type.to_csv(winners_with_type_fp, sep="\t", index=False)

    # Dropped genomes (not present in Cdb)
    input_genomes = set(lab["genome"].astype(str).tolist())
    used_genomes = set(df["genome"].astype(str).tolist())
    dropped = sorted(list(input_genomes - used_genomes))
    dropped_fp = os.path.join(report_dir, "dropped_by_drep.tsv")
    pd.DataFrame({"genome": dropped}).to_csv(dropped_fp, sep="\t", index=False)

    # ===== Split final dereplicated winner genomes into 3 folders =====
    # Classify winners by cluster_type; majority and winner_orig_label do not affect assignment.
    if not drep_rep_dir or not os.path.isdir(drep_rep_dir):
        die(f"Dereplicated genomes dir not found: {drep_rep_dir}")

    for sub in ["C1A", "C1B", "MIXED"]:
        os.makedirs(os.path.join(final_rep_dir, sub), exist_ok=True)

    # Use symlinks
    n_linked = 0
    missing = []

    for _, row in winners_with_type.iterrows():
        cluster_type = str(row["cluster_type"])
        winner = str(row["winner_genome"])

        # Map unexpected categories to MIXED so output remains limited to the three documented folders.
        if cluster_type not in {"C1A", "C1B", "MIXED"}:
            # Assign unexpected or missing types to MIXED to avoid dropping genomes.
            cluster_type = "MIXED"

        src = os.path.join(drep_rep_dir, winner)
        if not os.path.isfile(src):
            # Record winners that are absent from the dRep output directory.
            missing.append(winner)
            continue

        dst = os.path.join(final_rep_dir, cluster_type, winner)
        try:
            if os.path.islink(dst) or os.path.exists(dst):
                os.remove(dst)
            os.symlink(src, dst)
            n_linked += 1
        except Exception as e:
            die(f"Failed to link {src} -> {dst}: {e}")

    missing_fp = os.path.join(report_dir, "missing_winner_files.tsv")
    pd.DataFrame({"winner_genome_missing_in_dereplicated_genomes": missing}).to_csv(missing_fp, sep="\t", index=False)

    print("[DONE] Reports generated:")
    print("  -", cluster_summary_fp)
    print("  -", genome_cluster_map_fp)
    print("  -", mixed_cluster_genomes_fp)
    print("  -", winners_with_type_fp)
    print("  -", dropped_fp)
    print("  -", missing_fp)
    print("[DONE] Final reps split into 3 folders:")
    print("  -", os.path.join(final_rep_dir, "C1A"))
    print("  -", os.path.join(final_rep_dir, "C1B"))
    print("  -", os.path.join(final_rep_dir, "MIXED"))
    print(f"[DONE] Linked winner genomes: {n_linked}")

if __name__ == "__main__":
    main()
PYCODE

chmod +x "${PY}"

export CDB="${CDB}"
export WDB="${WDB}"
export LABELS_TSV="${LABELS_TSV}"
export REPORT_DIR="${REPORT_DIR}"
export DREP_REP_GENOMES_DIR="${DREP_OUT}/dereplicated_genomes"
export FINAL_REP_DIR="${FINAL_REP_DIR}"

echo "[INFO] Generating reports + split 3 folders..."
python3 "${PY}"

echo "============================================================"
echo "[ALL DONE]"
echo "Workdir:       ${WORKDIR}"
echo "dRep output:   ${DREP_OUT}"
echo "Reports:       ${REPORT_DIR}"
echo "Final reps:    ${FINAL_REP_DIR}"
echo ""
echo "Key outputs:"
echo "  - ${REPORT_DIR}/cluster_summary.tsv"
echo "  - ${REPORT_DIR}/genome_cluster_map.tsv"
echo "  - ${REPORT_DIR}/mixed_cluster_genomes.tsv"
echo "  - ${REPORT_DIR}/winners_with_cluster_type.tsv"
echo ""
echo "Final dereplicated winners split into:"
echo "  - ${FINAL_REP_DIR}/C1A"
echo "  - ${FINAL_REP_DIR}/C1B"
echo "  - ${FINAL_REP_DIR}/MIXED"
echo "============================================================"
