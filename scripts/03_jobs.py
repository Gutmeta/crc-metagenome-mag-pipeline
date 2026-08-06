import os
import pandas as pd

# 路径配置
sample_list_file = '/path/to/storage/data4/CRC_DATA/YachidaS_2019/SRR_Acc_List_CRC.txt'
output_job_dir = '3.DiTASiC_jobs/YachidaS_2019'
output_dir = '/path/to/storage/data3/CRC/YachidaS_2019/Results'
input_path = '/path/to/storage/data5/CRC_DATA/YachidaS_2019'
slurm_log_dir = '/path/to/storage/data3/CRC/YachidaS_2019/slurm_out'

# 创建脚本输出目录
os.makedirs(output_job_dir, exist_ok=True)
os.makedirs(slurm_log_dir, exist_ok=True)

# 读取样本列表
with open(sample_list_file, 'r') as f:
    samples = [line.strip() for line in f if line.strip()]

# 读取 seqkit_stats.txt 文件
seqkit_df = pd.read_csv("/path/to/storage/data4/CRC_DATA/YachidaS_2019/seqkit_stats.txt", sep=r'\s+', header=0)

# 提取样本ID（从 'file' 列中提取，样本ID为路径前部分）
seqkit_df['SampleID'] = seqkit_df['file'].str.split('/').str[0]

# 去掉 num_seqs 中的逗号并转换为整数
seqkit_df['num_seqs'] = seqkit_df['num_seqs'].str.replace(',', '').astype(int)


# 对每个样本的两个方向的 num_seqs 求和（双端加和）
total_reads_per_sample = seqkit_df.groupby('SampleID')['num_seqs'].sum()
min_reads = total_reads_per_sample.min()

# 遍历每个样本，生成对应的脚本
for sample in samples:
    script_path = os.path.join(output_job_dir, f"{sample}.sh")
    with open(script_path, 'w') as script:
        script.write(f"""#!/bin/bash
#SBATCH --job-name={sample}
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=64GB
#SBATCH --output={slurm_log_dir}/{sample}.out
#SBATCH --error={slurm_log_dir}/{sample}.err

ulimit -n 4096

export CHECKM_DATA_PATH="/path/to/storage/tools/CheckM"

# === Initialize Conda explicitly for non-interactive execution ===
CONDA_BIN="/path/to/conda/condabin/conda"
source /path/to/conda/etc/profile.d/conda.sh

# 让 conda 在非交互脚本中可用
eval "$("$CONDA_BIN" shell.bash hook)"

echo "Processing sample: {sample}"

fq_dir="{input_path}/{sample}/output/non_human_reads"
fq1=$(find "$fq_dir" -name "*_paired_1.fastq" | head -n 1)
fq2=$(find "$fq_dir" -name "*_paired_2.fastq" | head -n 1)

if [ ! -f "$fq1" ] || [ ! -f "$fq2" ]; then
    echo "Warning: Missing paired files for {sample}"
    exit 1
fi

mkdir -p "{output_dir}/DiTASiC/non_human_reads"
combined_fq="{output_dir}/DiTASiC/non_human_reads/{sample}_combined.fq"

mkdir -p "{output_dir}/DiTASiC/ditasic_mapping"

cat "$fq1" "$fq2" > "$combined_fq"

mkdir -p "{output_dir}/DiTASiC/ditasic_mapping_tmp"

set +u
conda activate kallisto046
set -u

ditasic_mapping.py -i {output_dir}/DiTASiC/kallisto_index \\
                   -l 100 \\
                   -t {output_dir}/DiTASiC/ditasic_mapping_tmp/tmp_{sample} \\
                   {output_dir}/DiTASiC/ref_paths.txt \\
                   "$combined_fq"

fq_base=$(basename "$combined_fq" .fq)
mapped_counts_file="${{fq_base}}_mapped_counts.npy"
total_counts_file="${{fq_base}}_total.npy"

if [ ! -f "$mapped_counts_file" ] || [ ! -f "$total_counts_file" ]; then
    echo "Error: ditasic_mapping failed for {sample}"
    exit 1
fi

mv "$mapped_counts_file" "$total_counts_file" "{output_dir}/DiTASiC/ditasic_mapping/"

mkdir -p "{output_dir}/DiTASiC/abundance"

rm "$combined_fq"

ditasic -r {output_dir}/DiTASiC/ref_paths.txt \\
        -a {output_dir}/DiTASiC/similarity_matrix.npy \\
        -x {output_dir}/DiTASiC/ditasic_mapping/$(basename "$mapped_counts_file") \\
        -n {output_dir}/DiTASiC/ditasic_mapping/$(basename "$total_counts_file") \\
        -o {output_dir}/DiTASiC/abundance/abundance_{sample}.txt

# 计算相对丰度
python - <<EOF
import pandas as pd

# 读取原始丰度文件
abundance_file = '{output_dir}/DiTASiC/abundance/abundance_{sample}.txt'
abundance_df = pd.read_csv(abundance_file, sep='\\t')

# 获取该样本的总reads数
total_reads = {total_reads_per_sample[sample]}

# 计算相对丰度
abundance_df['relative_abundance'] = abundance_df['count.estimate'] / total_reads

# 将相对丰度乘以最小reads数
abundance_df['adjusted_abundance'] = abundance_df['relative_abundance'] * {min_reads}

# 保存带有调整后的丰度的文件
adjusted_abundance_file = '{output_dir}/DiTASiC/abundance/adjusted_abundance_{sample}.txt'
abundance_df.to_csv(adjusted_abundance_file, sep='\\t', index=False)


EOF

""")

    print(f"✅ 脚本已生成: {script_path}")
