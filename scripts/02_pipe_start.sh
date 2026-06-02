#!/bin/bash

# 资源请求部分 (Slurm 指令)
#SBATCH --job-name=丰度计算CRC       # 作业名
#SBATCH --ntasks=1                  # 使用一个任务
#SBATCH --cpus-per-task=16           # 每个任务使用 8 个 CPU 核心
#SBATCH --mem=32GB                   # 请求 16GB 内存
#SBATCH --output=/path/to/storage/data4/CRC_DATA/YachidaS_2019/slurm_out/2.CRC.pipe.out
#SBATCH --error=/path/to/storage/data4/CRC_DATA/YachidaS_2019/slurm_out/2.CRC.pipe.err

bash ./2.pipe.sh /path/to/storage/data4/CRC_DATA/YachidaS_2019/SRR_Acc_List_CRC.txt /path/to/storage/data5/CRC_DATA/YachidaS_2019 /path/to/storage/data3/CRC/YachidaS_2019/Results
