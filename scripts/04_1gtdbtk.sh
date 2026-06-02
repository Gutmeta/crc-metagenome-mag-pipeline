#!/bin/bash
#SBATCH --job-name=gtdbtk_CRC
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=32GB
#SBATCH --output=/path/to/storage/data3/CRC/YachidaS_2019/slurm_out/gtdbtk.out
#SBATCH --error=/path/to/storage/data3/CRC/YachidaS_2019/slurm_out/gtdbtk.err


gtdbtk classify_wf  -x fa \
    --genome_dir /path/to/storage/data3/CRC/YachidaS_2019/Results/dRep_hq_bins_dir \
    --out_dir /path/to/crc-metagenome-mag-pipeline/classify_wf_out/YachidaS_2019 \
    --cpus 32

awk 'FNR==1 && NR!=1 {next} {print}' /path/to/crc-metagenome-mag-pipeline/classify_wf_out/YachidaS_2019/gtdbtk.*.summary.tsv > /path/to/crc-metagenome-mag-pipeline/classify_wf_out/YachidaS_2019/gtdbtk.summary.tsv