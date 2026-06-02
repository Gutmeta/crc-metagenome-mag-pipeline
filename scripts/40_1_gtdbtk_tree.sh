#!/bin/bash
#SBATCH --job-name=gtdbtk_tree
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=256
#SBATCH --mem=128GB
#SBATCH --output=/path/to/data2/CRC/CCDC2/slurm_out/gtdbtk_tree.out
#SBATCH --error=/path/to/data2/CRC/CCDC2/slurm_out/gtdbtk_tree.err

CONDA_BIN="/path/to/conda/condabin/conda"
eval "$("$CONDA_BIN" shell.bash hook)"
set +u
conda activate CRC
set -u

gtdbtk de_novo_wf \
  --genome_dir /path/to/data2/CRC/CCDC2/CC_TCG_genomes/ \
  --bacteria \
  --outgroup_taxon p__Patescibacteriota \
  --out_dir gtdbtk_denovo \
  -x fa \
  --cpus 256
