#!/bin/bash
#SBATCH --job-name=gtdbtk_FengQ_2015
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --mem=64GB
#SBATCH --nodelist=cn204
#SBATCH --chdir=/path/to/scratch/tmp_slurm
#SBATCH --output=/path/to/scratch/tmp_slurm/FengQ_2015/slurm_logs/gtdbtk_FengQ_2015.%j.out
#SBATCH --error=/path/to/scratch/tmp_slurm/FengQ_2015/slurm_logs/gtdbtk_FengQ_2015.%j.err

set -euo pipefail

CONDA_BIN="/path/to/conda/bin/conda"
eval "$($CONDA_BIN shell.bash hook)"
set +u
conda activate gtdbtk-2.4.1
set -u

export GTDBTK_DATA_PATH="/path/to/scratch/data3/gtdbtk_data/release226"

GENOME_DIR="/path/to/scratch/tmp_slurm/FengQ_2015/Results/dRep_hq_bins_dir"
OUT_DIR="/path/to/scratch/tmp_slurm/FengQ_2015/classify_wf_out/FengQ_2015"
SUMMARY_OUT="$OUT_DIR/gtdbtk.summary.tsv"

mkdir -p "$OUT_DIR"

gtdbtk classify_wf \
    --skip_ani_screen \
    -x fa \
    --genome_dir "$GENOME_DIR" \
    --out_dir "$OUT_DIR" \
    --cpus 128

awk 'FNR==1 && NR!=1 {next} {print}' "$OUT_DIR"/gtdbtk.*.summary.tsv > "$SUMMARY_OUT"
