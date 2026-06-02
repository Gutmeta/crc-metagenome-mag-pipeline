#!/bin/bash

ulimit -n 4096

export CHECKM_DATA_PATH="/path/to/databases/checkm_data"
# === Conda 初始化（用你的绝对路径，避免依赖 ~/.bashrc）===
CONDA_BIN="/path/to/conda/condabin/conda"
source /path/to/conda/etc/profile.d/conda.sh

# 让 conda 在非交互脚本中可用
eval "$("$CONDA_BIN" shell.bash hook)"

# 检查输入参数是否完整
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <sample_list> <output_dir>"
    exit 1
fi

# 获取输入参数
sample_list=$1
input_path=$2
output_dir=$3

mkdir -p ${output_dir}/dRep_hq_bins_dir/
if [ -d "${output_dir}/dRep_hq_bins_dir" ] && [ "$(ls -A ${output_dir}/dRep_hq_bins_dir)" ]; then
    rm -f ${output_dir}/dRep_hq_bins_dir/*
fi


shopt -s nullglob  # 启用空模式扩展
while read -r sample; do
    fa_files=(${input_path}/${sample}/output/dRep/dereplicated_genomes/*.fa)
    if [ ${#fa_files[@]} -gt 0 ]; then
        cp "${fa_files[@]}" ${output_dir}/dRep_hq_bins_dir/
    fi
done < ${sample_list}
shopt -u nullglob  # 关闭空模式扩展（可选）

# DiTASiC 计算丰度
set +u
conda activate CRC
set -u
# 生成 DiTASiC 相似性矩阵和索引
mkdir -p ${output_dir}/DiTASiC
find ${output_dir}/dRep_hq_bins_dir/ -name "*.fa" > ${output_dir}/DiTASiC/ref_paths.txt

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