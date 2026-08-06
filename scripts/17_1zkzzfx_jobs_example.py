import os

# 打开 sample_list.txt 文件，读取其中的 SRA 号
with open('/path/to/data1/FengQ_2015/SRR_Acc_List_CRC.txt', 'r') as f:
    samples = [line.strip() for line in f.readlines()]

# 创建 jobs 文件夹，如果不存在的话
os.makedirs('17.zkzzfx.jobs/FengQ_2015', exist_ok=True)

# 遍历每个 sample，生成对应的 Slurm 脚本
for sample in samples:
    # 定义 Slurm 脚本的内容
    slurm_script = f"""#!/bin/bash

# 资源请求部分 (Slurm 指令)
#SBATCH --job-name={sample}       # 作业名
#SBATCH --ntasks=1                  # 使用一个任务
#SBATCH --cpus-per-task=64           # 每个任务使用 64 个 CPU 核心
#SBATCH --mem=32GB                   # 请求 32GB 内存
#SBATCH --output=/path/to/data1/FengQ_2015/slurm_logs/{sample}.17.out
#SBATCH --error=/path/to/data1/FengQ_2015/slurm_logs/{sample}.17.err

bash /path/to/crc-metagenome-mag-pipeline/quality_control_assembly_binning.sh 64 \\
  "/path/to/data1/FengQ_2015/{sample}/{sample}_1.fastq.gz" \\
  "/path/to/data1/FengQ_2015/{sample}/{sample}_2.fastq.gz" \\
  "/path/to/data2/CRC/CCDC1_TCG/FengQ_2015/{sample}/output"
"""

    # 为每个 sample 创建一个对应的 slurm 脚本文件并保存在 jobs 文件夹中
    with open(f"17.zkzzfx.jobs/FengQ_2015/{sample}.sh", 'w') as out_file:
        out_file.write(slurm_script)

    print(f"Slurm 脚本 {sample}.sh 已生成并保存在 17.zkzzfx.jobs/FengQ_2015 文件夹中")
