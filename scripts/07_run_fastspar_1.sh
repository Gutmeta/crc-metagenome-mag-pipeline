#!/bin/bash
#SBATCH --job-name=fastspar_bootstrap_CRC_1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=128GB
#SBATCH --output=/path/to/storage/data3/CRC/YachidaS_2019/slurm_out/fastspar_bootstrap_CRC.out
#SBATCH --error=/path/to/storage/data3/CRC/YachidaS_2019/slurm_out/fastspar_bootstrap_CRC.err

set -e

# === Conda 初始化（用你的绝对路径，避免依赖 ~/.bashrc）===
CONDA_BIN="/path/to/conda/condabin/conda"
source ~/miniconda3/etc/profile.d/conda.sh

# 让 conda 在非交互脚本中可用
eval "$("$CONDA_BIN" shell.bash hook)"
set +u
conda activate fastspar-env
set -u

# 设置根目录，包含 CRC/control/adenoma 子目录
base_dir="/path/to/storage/data3/CRC/YachidaS_2019/Results/co_abundance_network_50"
groups=("CRC")
n_boot=1000

for group in "${groups[@]}"; do
  echo "🚀 处理组: $group"

  input_dir="${base_dir}/${group}"
  otu_table="${input_dir}/fastspar_input.tsv"
  corr_dir="${input_dir}/bootstrap_correlation"
  boot_dir="${input_dir}/bootstrap_counts"
  mkdir -p "$boot_dir" "$corr_dir"

  # 1️⃣ 原始相关性计算（用于和 bootstrap 分布比较）
  echo "📌 计算原始相关性矩阵"
  fastspar \
    --otu_table "$otu_table" \
    --correlation "${input_dir}/correlation.tsv" \
    --covariance "${input_dir}/covariance.tsv" \
    --iterations 1000 \
    --threads 64

  # 2️⃣ 生成 bootstrap OTU 表（1000 个）
  echo "📌 生成 bootstrap OTU 表"
  fastspar_bootstrap \
    --otu_table "$otu_table" \
    --number $n_boot \
    --prefix "${boot_dir}/otu"

  # 3️⃣ 批量计算 bootstrap 相关性矩阵
  echo "📌 批量计算 bootstrap 相关性"
  # 使用 GNU parallel 并发运行 fastspar
  parallel -j 64 fastspar \
    --otu_table {} \
    --correlation "${corr_dir}/cor_{/.}.tsv" \
    --covariance "${corr_dir}/cov_{/.}.tsv" \
    --iterations 1000 \
    --threads 1 ::: ${boot_dir}/otu_*.tsv

  # 4️⃣ 计算 p 值
  echo "📌 计算 p 值"
  fastspar_pvalues \
    --otu_table "$otu_table" \
    --correlation "${input_dir}/correlation.tsv" \
    --prefix "${corr_dir}/cor_otu_" \
    --permutations $n_boot \
    --outfile "${input_dir}/pvalues.tsv"

  echo "✅ $group 组完成"
done
