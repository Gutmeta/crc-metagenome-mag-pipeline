#!/bin/bash
#SBATCH --job-name=gtdbtk
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=128
#SBATCH --mem=64GB
#SBATCH --output=/path/to/data2/CRC/CCDC2/slurm_out/gtdbtk.out
#SBATCH --error=/path/to/data2/CRC/CCDC2/slurm_out/gtdbtk.err

CONDA_BIN="/path/to/conda/condabin/conda"
eval "$("$CONDA_BIN" shell.bash hook)"
set +u
conda activate CRC
set -u

gtdbtk classify_wf  -x fa \
    --genome_dir /path/to/data2/CRC/CCDC2/CC_TCG_genomes/ \
    --out_dir /path/to/crc-metagenome-mag-pipeline/classify_wf_out/CC_TCG \
    --cpus 128

awk 'FNR==1 && NR!=1 {next} {print}' /path/to/crc-metagenome-mag-pipeline/classify_wf_out/CC_TCG/gtdbtk.*.summary.tsv > /path/to/crc-metagenome-mag-pipeline/classify_wf_out/CC_TCG/gtdbtk.summary.tsv