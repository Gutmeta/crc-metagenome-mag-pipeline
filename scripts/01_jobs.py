import os

# 打开 sample_list.txt 文件，读取其中的 SRA 号
with open('sample_list.txt', 'r') as f:
    samples = [line.strip() for line in f.readlines()]

# 创建 jobs 文件夹，如果不存在的话
os.makedirs('1.jobs', exist_ok=True)

# 遍历每个 sample，生成对应的 Slurm 脚本
for sample in samples:
    # 定义 Slurm 脚本的内容
    slurm_script = f"""#!/bin/bash

# 资源请求部分 (Slurm 指令)
#SBATCH --job-name={sample}       # 作业名
#SBATCH --ntasks=1                  # 使用一个任务
#SBATCH --cpus-per-task=32           # 每个任务使用 32 个 CPU 核心
#SBATCH --mem=16GB                   # 请求 16GB 内存
#SBATCH --output=/path/to/storage/data/CRC/ThomasAM_2018/slurm_out/{sample}.out
#SBATCH --error=/path/to/storage/data/CRC/ThomasAM_2018/slurm_out/{sample}.err

bash /path/to/pipeline/quality_control_assembly_binning.sh 32 \\
  "/path/to/storage/data/CRC/ThomasAM_2018/{sample}/{sample}_1.fastq.gz" \\
  "/path/to/storage/data/CRC/ThomasAM_2018/{sample}/{sample}_2.fastq.gz" \\
  "/path/to/storage/data/CRC/ThomasAM_2018/{sample}/output"
"""

    # 为每个 sample 创建一个对应的 slurm 脚本文件并保存在 jobs 文件夹中
    with open(f"1.jobs/{sample}.sh", 'w') as out_file:
        out_file.write(slurm_script)

    print(f"Slurm 脚本 {sample}.sh 已生成并保存在 1.jobs 文件夹中")
