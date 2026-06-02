import pandas as pd
import os
import shutil
import re

# 文件路径设置
cluster_csv = "/path/to/storage/data3/CRC/YachidaS_2019/Results/co_abundance_network_75/C1_clusters_table.csv"
mapping_tsv = "/path/to/storage/data3/CRC/YachidaS_2019/Results/gtdbtk_id_mapping_with_unclassified.tsv"
source_dir = "/path/to/storage/data3/CRC/YachidaS_2019/Results/dRep_hq_bins_dir"
output_dir = "/path/to/storage/data3/CRC/CCDC1/HQMAG/YachidaS_2019_75"

# 创建输出目录（如果不存在）
os.makedirs(output_dir, exist_ok=True)

# 读取 cluster 表
cluster_df = pd.read_csv(cluster_csv)
# 提取 ID
# 提取最后一个空格后的字符串作为 ID
cluster_df['ID'] = cluster_df['Node'].apply(lambda x: x.strip().split()[-1])

# 读取映射表
mapping_df = pd.read_csv(mapping_tsv, sep="\t")

# 匹配 fa 文件名
merged_df = pd.merge(cluster_df, mapping_df, how='left', left_on='ID', right_on='ID')

# 检查哪些找到了匹配
not_found = merged_df[merged_df['fa_filename'].isnull()]
if not not_found.empty:
    print("⚠️ 以下ID未在映射文件中找到对应fa_filename：")
    print(not_found[['Node', 'ID']])

# 复制找到的 fasta 文件
for _, row in merged_df.dropna(subset=['fa_filename']).iterrows():
    fa_name = row['fa_filename']
    source_path = os.path.join(source_dir, fa_name)
    target_path = os.path.join(output_dir, fa_name)
    if os.path.exists(source_path):
        shutil.copy(source_path, target_path)
        print(f"✅ 已复制: {fa_name}")
    else:
        print(f"❌ 未找到文件: {source_path}")
