#!/bin/bash
#SBATCH --job-name=ANI_FengQ_2015
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=128GB
#SBATCH --output=/path/to/crc-metagenome-mag-pipeline/fastANI_results/FengQ_2015/slurm_fastANI.out
#SBATCH --error=/path/to/crc-metagenome-mag-pipeline/fastANI_results/FengQ_2015/slurm_fastANI.err

set -euo pipefail
shopt -s nullglob

C1B_DIR="/path/to/data2/CRC/YachidaS_2019/Results/C1_clusters/C1B"
CLUSTER1_DIR="/path/to/data2/CRC/CCDC1_TCG/FengQ_2015/Results/C1_clusters/cluster1"
CLUSTER2_DIR="/path/to/data2/CRC/CCDC1_TCG/FengQ_2015/Results/C1_clusters/cluster2"
OUTPUT_DIR="/path/to/crc-metagenome-mag-pipeline/fastANI_results/FengQ_2015"
FASTANI_BIN="/path/to/conda/envs/CRC/bin/fastANI"
THREADS="${SLURM_CPUS_PER_TASK:-64}"

mkdir -p "$OUTPUT_DIR"
: > "$OUTPUT_DIR/error_log.txt"
rm -f "$OUTPUT_DIR/cluster1_average_ani.txt" "$OUTPUT_DIR/cluster2_average_ani.txt"       "$OUTPUT_DIR/cluster1_label.txt" "$OUTPUT_DIR/cluster2_label.txt"

[ -x "$FASTANI_BIN" ] || { echo "[ERROR] fastANI 不存在: $FASTANI_BIN" >&2; exit 2; }
[ -d "$C1B_DIR" ] || { echo "[ERROR] C1B_DIR 不存在: $C1B_DIR" >&2; exit 2; }
[ -d "$CLUSTER1_DIR" ] || { echo "[ERROR] CLUSTER1_DIR 不存在: $CLUSTER1_DIR" >&2; exit 2; }
[ -d "$CLUSTER2_DIR" ] || { echo "[ERROR] CLUSTER2_DIR 不存在: $CLUSTER2_DIR" >&2; exit 2; }

C1B_LIST="$OUTPUT_DIR/c1b_genome_list.txt"
find "$C1B_DIR" -maxdepth 1 -type f -name '*.fa' | sort > "$C1B_LIST"
[ -s "$C1B_LIST" ] || { echo "[ERROR] C1B 参考列表为空: $C1B_LIST" >&2; exit 2; }

calculate_average_ani() {
    local cluster_dir="$1"
    local result_file="$2"
    local label="$3"
    local total_ani="0.0"
    local count=0
    local had_files=0

    for file in "$cluster_dir"/*.fa; do
        [ -f "$file" ] || continue
        had_files=1
        echo "[INFO] $label Processing $file"
        local tmp_out
        tmp_out=$(mktemp "$OUTPUT_DIR/${label}.XXXXXX.fastani.tsv")

        if ! "$FASTANI_BIN" -q "$file" --rl "$C1B_LIST" -o "$tmp_out" -t "$THREADS"; then
            echo "Error: fastANI failed for $file" >> "$OUTPUT_DIR/error_log.txt"
            rm -f "$tmp_out"
            continue
        fi

        local extracted_ani
        extracted_ani=$(python3 - "$tmp_out" <<'PYANI'
import sys
from pathlib import Path
p = Path(sys.argv[1])
vals = []
for line in p.read_text().splitlines():
    parts = line.split('	')
    if len(parts) >= 3:
        try:
            vals.append(float(parts[2]))
        except ValueError:
            pass
if vals:
    print(sum(vals) / len(vals))
PYANI
)
        rm -f "$tmp_out"

        if [[ -z "$extracted_ani" ]]; then
            echo "Error: Invalid ANI value for $file" >> "$OUTPUT_DIR/error_log.txt"
            continue
        fi

        total_ani=$(python3 - "$total_ani" "$extracted_ani" <<'PYSUM'
import sys
print(float(sys.argv[1]) + float(sys.argv[2]))
PYSUM
)
        count=$((count + 1))
    done

    if [[ "$had_files" -eq 0 ]]; then
        echo "Error: No .fa files found in $cluster_dir" >> "$OUTPUT_DIR/error_log.txt"
        exit 2
    fi

    if [[ "$count" -eq 0 ]]; then
        echo "Warning: No valid ANI values produced for $cluster_dir; fallback to 0.0" >> "$OUTPUT_DIR/error_log.txt"
        echo "0.0" > "$result_file"
        return 0
    fi

    python3 - "$total_ani" "$count" > "$result_file" <<'PYAVG'
import sys
print(f"{float(sys.argv[1]) / int(sys.argv[2]):.10f}")
PYAVG
}

cluster1_ani_file="$OUTPUT_DIR/cluster1_average_ani.txt"
cluster2_ani_file="$OUTPUT_DIR/cluster2_average_ani.txt"

calculate_average_ani "$CLUSTER1_DIR" "$cluster1_ani_file" cluster1
calculate_average_ani "$CLUSTER2_DIR" "$cluster2_ani_file" cluster2

cluster1_ani=$(cat "$cluster1_ani_file")
cluster2_ani=$(cat "$cluster2_ani_file")

if python3 - "$cluster1_ani" "$cluster2_ani" <<'PYCMP'
import sys
sys.exit(0 if float(sys.argv[1]) > float(sys.argv[2]) else 1)
PYCMP
then
    echo "Cluster 1 has higher ANI with C1B. Marking cluster1 as C1B and cluster2 as C1A."
    echo "cluster1 marked as C1B" > "$OUTPUT_DIR/cluster1_label.txt"
    echo "cluster2 marked as C1A" > "$OUTPUT_DIR/cluster2_label.txt"
else
    echo "Cluster 2 has higher ANI with C1B. Marking cluster2 as C1B and cluster1 as C1A."
    echo "cluster2 marked as C1B" > "$OUTPUT_DIR/cluster2_label.txt"
    echo "cluster1 marked as C1A" > "$OUTPUT_DIR/cluster1_label.txt"
fi

echo "Labeling complete. Results are saved in $OUTPUT_DIR."
