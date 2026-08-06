import os

# 路径配置
sample_list_file = '/path/to/storage/data2/FengQ_2015/SRR_Acc_List_CRC.txt'
output_job_dir   = '15_CCDC_1_DiTASiC_jobs/FengQ_2015'
output_dir       = '/path/to/storage/data3/CRC/CCDC1/FengQ_2015'
input_path       = '/path/to/storage/data2/FengQ_2015'
slurm_log_dir    = '/path/to/storage/data3/CRC/CCDC1/FengQ_2015/slurm_out'

# KneadData dependencies.
CONDA_BIN        = '/path/to/conda/condabin/conda'   # Explicit Conda executable.
CONDA_ENV        = 'kallisto046'
KNEADDATA_DB="/path/to/storage/tools/bowtie2-2.5.4-linux-x86_64/hg38_index/hg38_index"
TRIMMOMATIC_DIR="/path/to/conda/share/trimmomatic-0.39-2/"

# 创建脚本输出目录
os.makedirs(output_job_dir, exist_ok=True)
os.makedirs(slurm_log_dir, exist_ok=True)

# 读取样本列表
with open(sample_list_file, 'r') as f:
    samples = [line.strip() for line in f if line.strip()]

for sample in samples:
    script_path = os.path.join(output_job_dir, f"{sample}.sh")
    with open(script_path, 'w') as script:
        script.write(f"""#!/bin/bash
#SBATCH --job-name={sample}
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32GB
#SBATCH --output={slurm_log_dir}/{sample}.out
#SBATCH --error={slurm_log_dir}/{sample}.err

set -euo pipefail
ulimit -n 4096

export CHECKM_DATA_PATH="/path/to/checkm_data"

echo "=== Processing sample: {sample} ==="
threads="${{SLURM_CPUS_PER_TASK:-16}}"

# === Initialize Conda explicitly for non-interactive execution ===
set +u
CONDA_BIN="{CONDA_BIN}"
eval "$("$CONDA_BIN" shell.bash hook)"
set -u

# ---- 原始数据查找：在 input_path 根及下一层搜索，兼容 *_RmHost.1/.2, R1/R2, _1/_2, .1/.2 ----
input_root="{input_path}"

input_1=$(find "$input_root" -maxdepth 2 -type f \
  \\( -name "{sample}*_RmHost.1.fq.gz"  -o -name "{sample}*_RmHost.1.fastq.gz" \
     -name "{sample}*_RmHost.1.fq"     -o -name "{sample}*_RmHost.1.fastq"   \
     -o -name "{sample}*R1*.fastq.gz"  -o -name "{sample}*R1*.fq.gz"         \
     -o -name "{sample}*R1*.fastq"     -o -name "{sample}*R1*.fq"            \
     -o -name "{sample}*_1.fastq.gz"   -o -name "{sample}*_1.fq.gz"          \
     -o -name "{sample}*_1.fastq"      -o -name "{sample}*_1.fq"             \
     -o -name "{sample}*.1.fastq.gz"   -o -name "{sample}*.1.fq.gz"          \
     -o -name "{sample}*.1.fastq"      -o -name "{sample}*.1.fq" \\) \
  | head -n 1 || true)

input_2=$(find "$input_root" -maxdepth 2 -type f \
  \\( -name "{sample}*_RmHost.2.fq.gz"  -o -name "{sample}*_RmHost.2.fastq.gz" \
     -o -name "{sample}*_RmHost.2.fq"  -o -name "{sample}*_RmHost.2.fastq"   \
     -o -name "{sample}*R2*.fastq.gz"  -o -name "{sample}*R2*.fq.gz"         \
     -o -name "{sample}*R2*.fastq"     -o -name "{sample}*R2*.fq"            \
     -o -name "{sample}*_2.fastq.gz"   -o -name "{sample}*_2.fq.gz"          \
     -o -name "{sample}*_2.fastq"      -o -name "{sample}*_2.fq"             \
     -o -name "{sample}*.2.fastq.gz"   -o -name "{sample}*.2.fq.gz"          \
     -o -name "{sample}*.2.fastq"      -o -name "{sample}*.2.fq" \\) \
  | head -n 1 || true)

if [[ -z "$input_1" || -z "$input_2" ]]; then
  echo "[ERROR] Cannot find paired FASTQs for {sample} under {input_path} (searched depth<=2)" >&2
  exit 2
fi

echo "[INFO] Found inputs: $input_1 | $input_2"


# ---- 输出与前缀：直接用样本名，避免把 '.1' / '_RmHost' 带入 ----
base_name="{sample}"

# ==== kneaddata 质控+去宿主 ====
knead_out="{output_dir}/DiTASiC/non_human_reads"
mkdir -p "$knead_out"

KNEADDATA_DB="{KNEADDATA_DB}"
TRIMMOMATIC_DIR="{TRIMMOMATIC_DIR}"

# 已有结果则跳过（支持 .fastq 或 .fastq.gz）
paired1=$(find "$knead_out" -maxdepth 1 -type f \\( -name "${{base_name}}_paired_1.fastq" -o -name "${{base_name}}_paired_1.fastq.gz" \\) | head -n 1 || true)
paired2=$(find "$knead_out" -maxdepth 1 -type f \\( -name "${{base_name}}_paired_2.fastq" -o -name "${{base_name}}_paired_2.fastq.gz" \\) | head -n 1 || true)

if [ -z "$paired1" ] || [ -z "$paired2" ]; then
    echo "[kneaddata] Running for {sample} -> $knead_out/{sample}"
    extra_opts=()
    if [ -n "$TRIMMOMATIC_DIR" ]; then
        extra_opts+=(--trimmomatic "$TRIMMOMATIC_DIR")
    fi

    kneaddata \\
      --input1 "$input_1" \\
      --input2 "$input_2" \\
      --output "$knead_out" \\
      --output-prefix "$base_name" \\
      --reference-db "$KNEADDATA_DB" \\
      --threads "$threads" \\
      --remove-intermediate-output \\
      "${{extra_opts[@]}}"

    paired1=$(find "$knead_out" -maxdepth 1 -type f \\( -name "${{base_name}}_paired_1.fastq" -o -name "${{base_name}}_paired_1.fastq.gz" \\) | head -n 1 || true)
    paired2=$(find "$knead_out" -maxdepth 1 -type f \\( -name "${{base_name}}_paired_2.fastq" -o -name "${{base_name}}_paired_2.fastq.gz" \\) | head -n 1 || true)
fi

[ -s "$paired1" ] || {{ echo "[ERROR] KneadData 输出不存在：$paired1" >&2; exit 3; }}
[ -s "$paired2" ] || {{ echo "[ERROR] KneadData 输出不存在：$paired2" >&2; exit 3; }}
echo "[kneaddata] Done: $paired1 / $paired2"

# ===== DiTASiC 映射前：合并成单输入 =====
mkdir -p "{output_dir}/DiTASiC/ditasic_mapping" "{output_dir}/DiTASiC/ditasic_mapping_tmp" "{output_dir}/DiTASiC/abundance"
combined_fq="{output_dir}/DiTASiC/non_human_reads/{sample}_combined.fq"

# 把 paired_1 与 paired_2 依次写入同一个 fastq
if [[ "$paired1" == *.gz ]]; then zcat "$paired1" >  "$combined_fq"; else cat "$paired1" >  "$combined_fq"; fi
if [[ "$paired2" == *.gz ]]; then zcat "$paired2" >> "$combined_fq"; else cat "$paired2" >> "$combined_fq"; fi

# ===== DiTASiC 映射 =====
set +u
conda activate {CONDA_ENV}
set -u
ditasic_mapping.py -i {output_dir}/../DiTASiC_75/kallisto_index \\
                   -l 100 \\
                   -t {output_dir}/DiTASiC/ditasic_mapping_tmp/tmp_{sample} \\
                   {output_dir}/../DiTASiC_75/ref_paths.txt \\
                   "$combined_fq"

fq_base=$(basename "$combined_fq" .fq)
mapped_counts_file="${{fq_base}}_mapped_counts.npy"
total_counts_file="${{fq_base}}_total.npy"

if [ ! -f "$mapped_counts_file" ] || [ ! -f "$total_counts_file" ]; then
    echo "[ERROR] ditasic_mapping failed for {sample}" >&2
    exit 4
fi

mv "$mapped_counts_file" "$total_counts_file" "{output_dir}/DiTASiC/ditasic_mapping/"

# ===== DiTASiC 估丰度 =====
ditasic -r {output_dir}/../DiTASiC_75/ref_paths.txt \\
        -a {output_dir}/../DiTASiC_75/similarity_matrix.npy \\
        -x {output_dir}/DiTASiC/ditasic_mapping/$(basename "$mapped_counts_file") \\
        -n {output_dir}/DiTASiC/ditasic_mapping/$(basename "$total_counts_file") \\
        -o {output_dir}/DiTASiC/abundance/abundance_{sample}.txt
# ===== 成功后清理 non_human_reads 下以样本号开头的文件（仅文件）=====
cleanup_dir="{output_dir}/DiTASiC/non_human_reads"
echo "[CLEANUP] Removing files starting with '${{base_name}}' in $cleanup_dir"
if [ -d "$cleanup_dir" ]; then
  # 仅删除文件，保留目录；打印被删文件便于审计
  find "$cleanup_dir" -maxdepth 1 -type f -name "${{base_name}}*" -print -delete || true
fi

echo "[DONE] {sample} completed and cleaned up."
""")
    os.chmod(script_path, 0o755)
    print(f"✅ 脚本已生成: {script_path}")
