import os
import pandas as pd
import shutil

# 路径配置
fa_folder = '/path/to/data2/CRC/CCDC1/HQMAG/YachidaS_2019'
cluster_table_file = '/path/to/data2/CRC/YachidaS_2019/Results/co_abundance_network_75/C1_clusters_table.csv'
id_mapping_file = '/path/to/data2/CRC/YachidaS_2019/Results/gtdbtk_id_mapping_with_unclassified.tsv'
output_folder = '/path/to/data2/CRC/YachidaS_2019/Results/C1_clusters'

# 读取C1_clusters_table_2.csv文件
cluster_df = pd.read_csv(cluster_table_file)
cluster_df['ID'] = cluster_df['Node'].astype(str).str.split().str[-1]

# 读取gtdbtk_id_mapping_with_unclassified.tsv文件
id_mapping_df = pd.read_csv(id_mapping_file, sep='\t')

# 创建C1A和C1B文件夹
c1a_folder = os.path.join(output_folder, 'C1A')
c1b_folder = os.path.join(output_folder, 'C1B')

# 如果文件夹不存在则创建
os.makedirs(c1a_folder, exist_ok=True)
os.makedirs(c1b_folder, exist_ok=True)

# 合并两张表格，得到fa_filename与Group对应的信息
# 修改合并时的列名，ID列和Node列不同，需要使用 'ID' 来与 cluster_df 的 'Node' 合并
merged_df = pd.merge(
    id_mapping_df[['ID', 'fa_filename']],
    cluster_df[['ID', 'Group']],
    on='ID',
    how='inner'
)

# 遍历merged_df，按Group将.fa文件复制到对应的文件夹
for _, row in merged_df.iterrows():
    fa_filename = row['fa_filename']
    group = row['Group']
    
    # 构造文件路径
    fa_file_path = os.path.join(fa_folder, fa_filename)
    
    # 确保.fa文件存在
    if os.path.exists(fa_file_path):
        # 根据Group将文件复制到相应的文件夹
        if group == 'C1A':
            shutil.copy(fa_file_path, c1a_folder)
        elif group == 'C1B':
            shutil.copy(fa_file_path, c1b_folder)
    else:
        print(f"Warning: {fa_file_path} not found.")

print("文件整理完成。")
