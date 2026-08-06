#!/usr/bin/env python3
import argparse
import csv
import os
from dataclasses import dataclass
from pathlib import Path


DATASET_NAME = "WirbelJ_2019"

# Input metadata
run_group_file = Path("/path/to/data3/CRC_DATA/CCDC2_data/WirbelJ_2019/PRJEB27928_PE_runs_grouped.txt")
group_file = Path("WirbelJ_2019_CRC_Group.txt")

# Paths
output_job_dir = Path("41_CC_TCG_DiTASiC_jobs/WirbelJ_2019_optimized")
output_dir = Path("/path/to/data2/CRC/CCDC2_val/WirbelJ_2019")
input_path = Path("/path/to/data3/CRC_DATA/CCDC2_data/WirbelJ_2019")
slurm_log_dir = output_dir / "slurm_out_optimized"
scratch_root = Path("/tmp/WirbelJ_2019_DiTASiC")
discard_samples_file = Path("WirbelJ_2019_DiTASiC_discard_low_assignment.txt")

# kneaddata / DiTASiC dependencies
CONDA_BIN = "/path/to/conda/condabin/conda"
CONDA_ENV = "CRC"
CRC_ENV_PREFIX = "/path/to/conda/envs/CRC"
KNEADDATA_DB = "/path/to/databases/hg38/index/hg38_index"
TRIMMOMATIC_DIR = "/path/to/conda/envs/CRC/share/trimmomatic-0.40-0"
REFERENCE_ROOT = str(output_dir / "../DiTASiC")
DITASIC_BIN = "/path/to/ditasic/ditasic"
DITASIC_MAPPING_PY = "/path/to/ditasic/ditasic_mapping.py"
KALLISTO_BIN_DIR = "/path/to/ditasic/bin"

DEFAULT_THREADS = 32
DEFAULT_MEM = "64G"
DEFAULT_KNEADDATA_MAX_MEMORY = "30000M"
DEFAULT_ARRAY_CONCURRENCY = 4


@dataclass(frozen=True)
class SampleInfo:
    sample: str
    group: str
    runs: tuple[str, ...]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Generate DiTASiC scripts for WirbelJ_2019. "
            "Each PE grouped row is treated as one biological sample."
        )
    )
    parser.add_argument(
        "--node",
        help="Optional Slurm node constraint; omitted by default.",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=DEFAULT_THREADS,
        help="CPUs per task and tool threads. Default: 32",
    )
    parser.add_argument("--mem", default=DEFAULT_MEM, help="Slurm memory per array task. Default: 64G")
    parser.add_argument(
        "--kneaddata-max-memory",
        default=DEFAULT_KNEADDATA_MAX_MEMORY,
        help="Passed to kneaddata --max-memory. Default: 30000M",
    )
    parser.add_argument(
        "--array-concurrency",
        type=int,
        default=DEFAULT_ARRAY_CONCURRENCY,
        help="Slurm array concurrency. Default: 4.",
    )
    parser.add_argument(
        "--discard-samples-file",
        type=Path,
        default=discard_samples_file,
        help="Samples to skip because DiTASiC assigned ratio is below threshold.",
    )
    return parser.parse_args()


def split_runs(value: str) -> list[str]:
    return [item.strip() for item in str(value).split(";") if item.strip()]


def fastq_pair_exists(run: str) -> bool:
    return (input_path / run / f"{run}_1.fastq.gz").is_file() and (
        input_path / run / f"{run}_2.fastq.gz"
    ).is_file()


def load_samples(path: Path) -> list[str]:
    if not path.exists():
        return []
    return [
        line.strip()
        for line in path.read_text().splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    ]


def load_group_by_run(path: Path) -> dict[str, str]:
    if not path.exists():
        raise FileNotFoundError(f"Group file not found: {path}")

    group_by_run: dict[str, str] = {}
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or "Sample" not in reader.fieldnames or "Group" not in reader.fieldnames:
            raise ValueError(f"{path} must contain Sample and Group columns")
        for row in reader:
            group = row["Group"].strip()
            if not group:
                continue
            for run in split_runs(row["Sample"]):
                previous = group_by_run.get(run)
                if previous is not None and previous != group:
                    raise ValueError(f"Conflicting group labels for {run}: {previous} vs {group}")
                group_by_run[run] = group
    return group_by_run


def load_sample_infos() -> tuple[list[SampleInfo], int, int]:
    if not run_group_file.exists():
        raise FileNotFoundError(f"Run group file not found: {run_group_file}")

    group_by_run = load_group_by_run(group_file)
    samples: list[SampleInfo] = []
    ignored_missing_pe_runs = 0
    total_downloaded_runs = 0

    for line_no, line in enumerate(run_group_file.read_text().splitlines(), start=1):
        pe_runs = split_runs(line)
        if not pe_runs:
            continue

        downloaded_runs = tuple(run for run in pe_runs if fastq_pair_exists(run))
        ignored_missing_pe_runs += len(pe_runs) - len(downloaded_runs)
        total_downloaded_runs += len(downloaded_runs)
        if not downloaded_runs:
            continue

        labels = {group_by_run[run] for run in downloaded_runs if run in group_by_run}
        missing_labels = [run for run in downloaded_runs if run not in group_by_run]
        if missing_labels:
            raise ValueError(
                f"Missing group label for downloaded PE runs on line {line_no}: "
                + ";".join(missing_labels)
            )
        if len(labels) != 1:
            raise ValueError(
                f"Downloaded PE runs on line {line_no} have inconsistent labels: "
                + ";".join(downloaded_runs)
            )

        sample = downloaded_runs[0]
        samples.append(SampleInfo(sample=sample, group=labels.pop(), runs=downloaded_runs))

    seen_samples: set[str] = set()
    duplicate_samples = [info.sample for info in samples if info.sample in seen_samples or seen_samples.add(info.sample)]
    if duplicate_samples:
        raise ValueError("Duplicate sample IDs after PE filtering: " + ", ".join(duplicate_samples))

    return samples, total_downloaded_runs, ignored_missing_pe_runs


def existing_abundance_samples(abundance_dir: Path) -> set[str]:
    return {
        path.name.removeprefix("abundance_").removesuffix(".txt")
        for path in abundance_dir.glob("abundance_*.txt")
        if path.stat().st_size > 0
    }


def write_executable(path: Path, text: str) -> None:
    path.write_text(text)
    os.chmod(path, 0o755)


def replace_tokens(template: str, values: dict[str, object]) -> str:
    text = template
    for key, value in values.items():
        text = text.replace(f"__{key}__", str(value))
    return text


def shell_run_list(runs: tuple[str, ...]) -> str:
    return ";".join(runs)


SAMPLE_SCRIPT_TEMPLATE = r"""#!/usr/bin/env bash
#SBATCH --job-name=__SAMPLE__
__NODE_DIRECTIVE__
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=__THREADS__
#SBATCH --mem=__MEM__
#SBATCH --output=__SLURM_LOG_DIR__/__SAMPLE__.%j.out
#SBATCH --error=__SLURM_LOG_DIR__/__SAMPLE__.%j.err

set -Eeuo pipefail
trap 'echo "[ERROR] command failed at line $LINENO: $BASH_COMMAND" >&2' ERR
umask 002
ulimit -n 4096

SAMPLE="__SAMPLE__"
GROUP="__GROUP__"
RUN_LIST="__RUN_LIST__"
THREADS="${SLURM_CPUS_PER_TASK:-__THREADS__}"
KNEADDATA_MAX_MEMORY="${KNEADDATA_MAX_MEMORY:-__KNEADDATA_MAX_MEMORY__}"
INPUT_ROOT="__INPUT_PATH__"
OUTPUT_DIR="__OUTPUT_DIR__"
REFERENCE_ROOT="__REFERENCE_ROOT__"
SCRATCH_ROOT="${SCRATCH_ROOT:-__SCRATCH_ROOT__}"
SCRATCH_DIR="${SCRATCH_ROOT}/${SAMPLE}_${SLURM_JOB_ID:-manual_$$}_${SLURM_ARRAY_TASK_ID:-0}"
SLURM_LOG_DIR="__SLURM_LOG_DIR__"

CONDA_BIN="__CONDA_BIN__"
CONDA_ENV="__CONDA_ENV__"
CRC_ENV_PREFIX="__CRC_ENV_PREFIX__"
KNEADDATA_DB="__KNEADDATA_DB__"
TRIMMOMATIC_DIR="__TRIMMOMATIC_DIR__"
DITASIC_BIN="__DITASIC_BIN__"
DITASIC_MAPPING_PY="__DITASIC_MAPPING_PY__"
KALLISTO_BIN_DIR="__KALLISTO_BIN_DIR__"

RAW1="${SCRATCH_DIR}/${SAMPLE}_raw_1.fastq.gz"
RAW2="${SCRATCH_DIR}/${SAMPLE}_raw_2.fastq.gz"
KNEAD_DIR="${SCRATCH_DIR}/non_human_reads"
MAPPING_TMP_DIR="${SCRATCH_DIR}/ditasic_mapping_tmp"
COMBINED_FQ="${SCRATCH_DIR}/${SAMPLE}_combined.fq"
LOCAL_MAPPED="${SCRATCH_DIR}/${SAMPLE}_combined_mapped_counts.npy"
LOCAL_TOTAL="${SCRATCH_DIR}/${SAMPLE}_combined_total.npy"
FINAL_MAPPING_DIR="${OUTPUT_DIR}/DiTASiC/ditasic_mapping"
FINAL_ABUNDANCE_DIR="${OUTPUT_DIR}/DiTASiC/abundance"
FINAL_MAPPED="${FINAL_MAPPING_DIR}/${SAMPLE}_combined_mapped_counts.npy"
FINAL_TOTAL="${FINAL_MAPPING_DIR}/${SAMPLE}_combined_total.npy"
ABUNDANCE_FILE="${FINAL_ABUNDANCE_DIR}/abundance_${SAMPLE}.txt"

cleanup() {
  status=$?
  if [[ -n "${KNEAD_DIR:-}" && -f "${KNEAD_DIR}/${SAMPLE}.log" ]]; then
    mkdir -p "${SLURM_LOG_DIR}" || true
    cp -f "${KNEAD_DIR}/${SAMPLE}.log" "${SLURM_LOG_DIR}/${SAMPLE}.${SLURM_JOB_ID:-manual}.kneaddata.log" || true
  fi
  if [[ -n "${SCRATCH_DIR:-}" && -d "${SCRATCH_DIR}" && "${SCRATCH_DIR}" == "${SCRATCH_ROOT}/${SAMPLE}"_* ]]; then
    echo "[CLEANUP] removing scratch: ${SCRATCH_DIR}"
    rm -rf "${SCRATCH_DIR}"
  fi
  exit "${status}"
}
trap cleanup EXIT

echo "[INFO] start: $(date '+%F %T')"
echo "[INFO] node=$(hostname) sample=${SAMPLE} group=${GROUP} threads=${THREADS} kneaddata_max_memory=${KNEADDATA_MAX_MEMORY}"
echo "[INFO] runs=${RUN_LIST}"
echo "[INFO] scratch=${SCRATCH_DIR}"
echo "[INFO] abundance=${ABUNDANCE_FILE}"

if [[ -s "${ABUNDANCE_FILE}" ]]; then
  echo "[SKIP] abundance already exists: ${ABUNDANCE_FILE}"
  exit 0
fi

mkdir -p "${KNEAD_DIR}" "${MAPPING_TMP_DIR}" "${FINAL_MAPPING_DIR}" "${FINAL_ABUNDANCE_DIR}"

IFS=';' read -r -a RUNS <<< "${RUN_LIST}"
if [[ "${#RUNS[@]}" -eq 0 ]]; then
  echo "[ERROR] no downloaded PE runs for ${SAMPLE}" >&2
  exit 2
fi

: > "${RAW1}"
: > "${RAW2}"
for RUN in "${RUNS[@]}"; do
  INPUT1="${INPUT_ROOT}/${RUN}/${RUN}_1.fastq.gz"
  INPUT2="${INPUT_ROOT}/${RUN}/${RUN}_2.fastq.gz"
  [[ -s "${INPUT1}" ]] || { echo "[ERROR] missing input1 for ${RUN}: ${INPUT1}" >&2; exit 2; }
  [[ -s "${INPUT2}" ]] || { echo "[ERROR] missing input2 for ${RUN}: ${INPUT2}" >&2; exit 2; }
  echo "[INFO] adding run ${RUN}"
  cat "${INPUT1}" >> "${RAW1}"
  cat "${INPUT2}" >> "${RAW2}"
done

[[ -s "${RAW1}" ]] || { echo "[ERROR] raw R1 is empty: ${RAW1}" >&2; exit 2; }
[[ -s "${RAW2}" ]] || { echo "[ERROR] raw R2 is empty: ${RAW2}" >&2; exit 2; }
[[ -x "${CONDA_BIN}" ]] || { echo "[ERROR] missing conda: ${CONDA_BIN}" >&2; exit 3; }
[[ -d "${TRIMMOMATIC_DIR}" ]] || { echo "[ERROR] missing trimmomatic dir: ${TRIMMOMATIC_DIR}" >&2; exit 4; }
([[ -e "${KNEADDATA_DB}.1.bt2" ]] || [[ -e "${KNEADDATA_DB}.1.bt2l" ]]) || { echo "[ERROR] missing kneaddata DB prefix: ${KNEADDATA_DB}" >&2; exit 5; }
[[ -r "${REFERENCE_ROOT}/kallisto_index" ]] || { echo "[ERROR] missing kallisto_index: ${REFERENCE_ROOT}/kallisto_index" >&2; exit 6; }
[[ -r "${REFERENCE_ROOT}/ref_paths.txt" ]] || { echo "[ERROR] missing ref_paths.txt: ${REFERENCE_ROOT}/ref_paths.txt" >&2; exit 7; }
[[ -r "${REFERENCE_ROOT}/similarity_matrix.npy" ]] || { echo "[ERROR] missing similarity_matrix.npy: ${REFERENCE_ROOT}/similarity_matrix.npy" >&2; exit 8; }
[[ -r "${DITASIC_MAPPING_PY}" ]] || { echo "[ERROR] missing ditasic_mapping.py: ${DITASIC_MAPPING_PY}" >&2; exit 9; }
[[ -x "${DITASIC_BIN}" ]] || { echo "[ERROR] missing ditasic: ${DITASIC_BIN}" >&2; exit 10; }

eval "$("${CONDA_BIN}" shell.bash hook)"
set +u
conda activate "${CONDA_ENV}"
set -u
export PATH="${CRC_ENV_PREFIX}/bin:${KALLISTO_BIN_DIR}:/path/to/ditasic:${PATH}"
export CHECKM_DATA_PATH="/path/to/checkm_data"

command -v kneaddata >/dev/null || { echo "[ERROR] kneaddata not in PATH" >&2; exit 11; }
command -v kallisto >/dev/null || { echo "[ERROR] kallisto not in PATH" >&2; exit 12; }

rm -f "${FINAL_MAPPED}" "${FINAL_TOTAL}"

echo "[TIME] kneaddata start: $(date '+%F %T')"
kneaddata \
  --input1 "${RAW1}" \
  --input2 "${RAW2}" \
  --output "${KNEAD_DIR}" \
  --output-prefix "${SAMPLE}" \
  --reference-db "${KNEADDATA_DB}" \
  --threads "${THREADS}" \
  --max-memory "${KNEADDATA_MAX_MEMORY}" \
  --remove-intermediate-output \
  --trimmomatic "${TRIMMOMATIC_DIR}"
echo "[TIME] kneaddata end: $(date '+%F %T')"

paired1=$(find "${KNEAD_DIR}" -maxdepth 1 -type f \( -name "${SAMPLE}_paired_1.fastq" -o -name "${SAMPLE}_paired_1.fastq.gz" \) | head -n 1 || true)
paired2=$(find "${KNEAD_DIR}" -maxdepth 1 -type f \( -name "${SAMPLE}_paired_2.fastq" -o -name "${SAMPLE}_paired_2.fastq.gz" \) | head -n 1 || true)

[[ -s "${paired1}" ]] || { echo "[ERROR] missing kneaddata output paired1: ${paired1}" >&2; exit 13; }
[[ -s "${paired2}" ]] || { echo "[ERROR] missing kneaddata output paired2: ${paired2}" >&2; exit 14; }

echo "[TIME] combine start: $(date '+%F %T')"
if [[ "${paired1}" == *.gz ]]; then zcat "${paired1}" > "${COMBINED_FQ}"; else cat "${paired1}" > "${COMBINED_FQ}"; fi
if [[ "${paired2}" == *.gz ]]; then zcat "${paired2}" >> "${COMBINED_FQ}"; else cat "${paired2}" >> "${COMBINED_FQ}"; fi
[[ -s "${COMBINED_FQ}" ]] || { echo "[ERROR] combined FASTQ is empty: ${COMBINED_FQ}" >&2; exit 15; }
echo "[TIME] combine end: $(date '+%F %T')"

echo "[TIME] ditasic_mapping start: $(date '+%F %T')"
pushd "${SCRATCH_DIR}" >/dev/null
python "${DITASIC_MAPPING_PY}" \
  -i "${REFERENCE_ROOT}/kallisto_index" \
  -l 100 \
  -t "${MAPPING_TMP_DIR}/tmp_${SAMPLE}" \
  "${REFERENCE_ROOT}/ref_paths.txt" \
  "${COMBINED_FQ}"
popd >/dev/null
echo "[TIME] ditasic_mapping end: $(date '+%F %T')"

[[ -f "${LOCAL_MAPPED}" ]] || { echo "[ERROR] missing mapped counts: ${LOCAL_MAPPED}" >&2; exit 16; }
[[ -f "${LOCAL_TOTAL}" ]] || { echo "[ERROR] missing total counts: ${LOCAL_TOTAL}" >&2; exit 17; }
mv -f "${LOCAL_MAPPED}" "${FINAL_MAPPED}"
mv -f "${LOCAL_TOTAL}" "${FINAL_TOTAL}"

echo "[TIME] ditasic start: $(date '+%F %T')"
"${DITASIC_BIN}" \
  -r "${REFERENCE_ROOT}/ref_paths.txt" \
  -a "${REFERENCE_ROOT}/similarity_matrix.npy" \
  -x "${FINAL_MAPPED}" \
  -n "${FINAL_TOTAL}" \
  -o "${ABUNDANCE_FILE}"
echo "[TIME] ditasic end: $(date '+%F %T')"

[[ -s "${ABUNDANCE_FILE}" ]] || { echo "[ERROR] missing abundance output: ${ABUNDANCE_FILE}" >&2; exit 18; }

echo "[DONE] ${SAMPLE} completed: $(date '+%F %T')"
"""


ARRAY_SCRIPT_TEMPLATE = r"""#!/usr/bin/env bash
#SBATCH --job-name=WirbelJ19_DiTASiC
__NODE_DIRECTIVE__
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=__THREADS__
#SBATCH --mem=__MEM__
#SBATCH --array=0-__LAST_INDEX__%__ARRAY_CONCURRENCY__
#SBATCH --output=__SLURM_LOG_DIR__/array_%A_%a.out
#SBATCH --error=__SLURM_LOG_DIR__/array_%A_%a.err

set -Eeuo pipefail
trap 'echo "[ERROR] array task failed at line $LINENO: $BASH_COMMAND" >&2' ERR

SAMPLE_LIST="__SAMPLE_LIST__"
JOB_DIR="__JOB_DIR__"
export SCRATCH_ROOT="${SCRATCH_ROOT:-__SCRATCH_ROOT__}"
export KNEADDATA_MAX_MEMORY="${KNEADDATA_MAX_MEMORY:-__KNEADDATA_MAX_MEMORY__}"

task_id="${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID is required}"
sample=$(sed -n "$((task_id + 1))p" "${SAMPLE_LIST}")

if [[ -z "${sample}" ]]; then
  echo "[ERROR] no sample for array task ${task_id}" >&2
  exit 2
fi

script="${JOB_DIR}/${sample}.sh"
[[ -x "${script}" ]] || { echo "[ERROR] missing executable sample script: ${script}" >&2; exit 3; }

launch_dir="${SCRATCH_ROOT}/array_launch_scripts/${SLURM_ARRAY_JOB_ID:-manual}_${task_id}_$$"
launch_script="${launch_dir}/${sample}.sh"

cleanup_array() {
  status=$?
  if [[ -n "${launch_dir:-}" && -d "${launch_dir}" ]]; then
    rm -rf "${launch_dir}"
  fi
  exit "${status}"
}
trap cleanup_array EXIT

mkdir -p "${launch_dir}"
cp -f "${script}" "${launch_script}"
chmod +x "${launch_script}"

echo "[INFO] array_job=${SLURM_ARRAY_JOB_ID:-NA} task=${task_id} sample=${sample}"
echo "[INFO] script=${script}"
echo "[INFO] launch_script=${launch_script}"
echo "[INFO] cpus=${SLURM_CPUS_PER_TASK:-NA} scratch_root=${SCRATCH_ROOT} kneaddata_max_memory=${KNEADDATA_MAX_MEMORY}"

bash "${launch_script}"
"""


def sample_script_text(info: SampleInfo, args: argparse.Namespace) -> str:
    return replace_tokens(
        SAMPLE_SCRIPT_TEMPLATE,
        {
            "SAMPLE": info.sample,
            "GROUP": info.group,
            "RUN_LIST": shell_run_list(info.runs),
            "NODE_DIRECTIVE": f"#SBATCH --nodelist={args.node}" if args.node else "",
            "THREADS": args.threads,
            "MEM": args.mem,
            "SLURM_LOG_DIR": slurm_log_dir,
            "KNEADDATA_MAX_MEMORY": args.kneaddata_max_memory,
            "INPUT_PATH": input_path,
            "OUTPUT_DIR": output_dir,
            "REFERENCE_ROOT": REFERENCE_ROOT,
            "SCRATCH_ROOT": scratch_root,
            "CONDA_BIN": CONDA_BIN,
            "CONDA_ENV": CONDA_ENV,
            "CRC_ENV_PREFIX": CRC_ENV_PREFIX,
            "KNEADDATA_DB": KNEADDATA_DB,
            "TRIMMOMATIC_DIR": TRIMMOMATIC_DIR,
            "DITASIC_BIN": DITASIC_BIN,
            "DITASIC_MAPPING_PY": DITASIC_MAPPING_PY,
            "KALLISTO_BIN_DIR": KALLISTO_BIN_DIR,
        },
    )


def array_script_text(missing_samples_file: Path, args: argparse.Namespace, last_index: int) -> str:
    return replace_tokens(
        ARRAY_SCRIPT_TEMPLATE,
        {
            "NODE_DIRECTIVE": f"#SBATCH --nodelist={args.node}" if args.node else "",
            "THREADS": args.threads,
            "MEM": args.mem,
            "LAST_INDEX": last_index,
            "ARRAY_CONCURRENCY": args.array_concurrency,
            "SLURM_LOG_DIR": slurm_log_dir,
            "SAMPLE_LIST": missing_samples_file.resolve(),
            "JOB_DIR": output_job_dir.resolve(),
            "SCRATCH_ROOT": scratch_root,
            "KNEADDATA_MAX_MEMORY": args.kneaddata_max_memory,
        },
    )


def write_manifest(path: Path, samples: list[SampleInfo]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(["Sample", "Group", "Runs", "NumRuns"])
        for info in samples:
            writer.writerow([info.sample, info.group, ";".join(info.runs), len(info.runs)])


def main() -> None:
    args = parse_args()
    output_job_dir.mkdir(parents=True, exist_ok=True)
    slurm_log_dir.mkdir(parents=True, exist_ok=True)

    sample_infos, total_downloaded_runs, ignored_missing_pe_runs = load_sample_infos()
    manifest_file = output_job_dir / "WirbelJ_2019_DiTASiC_sample_manifest.tsv"
    write_manifest(manifest_file, sample_infos)

    abundance_dir = output_dir / "DiTASiC" / "abundance"
    done = existing_abundance_samples(abundance_dir)
    discarded = set(load_samples(args.discard_samples_file))
    missing_infos = [
        info for info in sample_infos if info.sample not in done and info.sample not in discarded
    ]

    missing_samples_file = output_job_dir / "WirbelJ_2019_DiTASiC_samples.txt"
    missing_samples_file.write_text("".join(f"{info.sample}\n" for info in missing_infos))

    for info in missing_infos:
        script_path = output_job_dir / f"{info.sample}.sh"
        write_executable(script_path, sample_script_text(info, args))

    if missing_infos:
        array_script = output_job_dir / "submit_array_WirbelJ_2019_not_done.sh"
        write_executable(array_script, array_script_text(missing_samples_file, args, len(missing_infos) - 1))
    else:
        array_script = None

    print(f"dataset={DATASET_NAME}")
    print(f"total_sample_groups={len(sample_infos)}")
    print(f"downloaded_pe_runs_used={total_downloaded_runs}")
    print(f"ignored_missing_pe_runs={ignored_missing_pe_runs}")
    print(f"existing_abundance={len(done)}")
    print(f"discarded_low_assignment={len(discarded)}")
    print(f"missing_abundance={len(missing_infos)}")
    print(f"output_job_dir={output_job_dir}")
    print(f"manifest_file={manifest_file}")
    print(f"missing_samples_file={missing_samples_file}")
    if array_script:
        print(f"array_script={array_script}")
        print(f"submit_command=sbatch {array_script}")
    else:
        print("array_script=NONE (all samples already have abundance outputs)")


if __name__ == "__main__":
    main()
