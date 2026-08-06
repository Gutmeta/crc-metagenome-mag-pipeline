#!/bin/bash

# 资源请求部分 (Slurm 指令)
#SBATCH --job-name=丰度计算CRC       # 作业名
#SBATCH --ntasks=1                  # 使用一个任务
#SBATCH --cpus-per-task=64           # 每个任务使用 64 个 CPU 核心
#SBATCH --mem=128GB                   # 请求 128GB 内存
#SBATCH --output=/path/to/data2/CRC/CCDC1/slurm_out/14.pipe.out
#SBATCH --error=/path/to/data2/CRC/CCDC1/slurm_out/14.pipe.err

# 获取输入参数
output_dir=/path/to/data2/CRC/CCDC1/
ulimit -n 4096

# === Initialize Conda explicitly for non-interactive execution ===
CONDA_BIN="/path/to/conda/bin/conda"
source /path/to/conda/etc/profile.d/conda.sh

# 让 conda 在非交互脚本中可用
eval "$("$CONDA_BIN" shell.bash hook)"

set +u
conda activate CRC
set -u
# DiTASiC 计算丰度
#生成 DiTASiC 相似性矩阵和索引
mkdir -p ${output_dir}/DiTASiC
find ${output_dir}/HQMAG/YachidaS_2019/ -name "*.fa" > ${output_dir}/DiTASiC/ref_paths.txt

# 手动生成 kallisto 索引（关键步骤）
echo "Generating kallisto index manually..."
kallisto index -i ${output_dir}/DiTASiC/kallisto_index \
               --make-unique \
               $(cat ${output_dir}/DiTASiC/ref_paths.txt)

# 生成相似性矩阵（指定手动生成的索引）
ditasic_matrix.py -l 100 \
               ${output_dir}/DiTASiC/ref_paths.txt \
               -t ${output_dir}/DiTASiC/temp \
               -i ${output_dir}/DiTASiC/kallisto_index \
               -o ${output_dir}/DiTASiC/similarity_matrix.npy

# 检查是否生成矩阵和索引
if [ ! -f "${output_dir}/DiTASiC/similarity_matrix.npy" ] || [ ! -f "${output_dir}/DiTASiC/kallisto_index" ]; then
    echo "Error: ditasic_matrix failed to generate similarity matrix or index."
    exit 1
fi
