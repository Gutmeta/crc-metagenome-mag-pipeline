import os
import pandas as pd
from glob import glob

# ==== 配置路径 ====
abundance_dir = "/path/to/data2/CRC/CCDC1_TCG/FengQ_2015/Results/DiTASiC/abundance"
gtdb_summary_file = "/path/to/crc-metagenome-mag-pipeline/classify_wf_out/FengQ_2015/gtdbtk.summary.tsv"
output_file = "/path/to/data2/CRC/CCDC1_TCG/FengQ_2015/Results/merged_abundance_with_taxonomy.tsv"

# ==== 读取 GTDB-Tk summary 文件 ====
taxonomy_df = pd.read_csv(gtdb_summary_file, sep="\t", usecols=["user_genome", "classification"])
taxonomy_df["taxa.name"] = taxonomy_df["user_genome"].astype(str) + ".fa"
taxonomy_df = taxonomy_df[["taxa.name", "classification"]]
taxonomy_df.columns = ["taxa.name", "taxonomy"]

# ==== 合并 abundance_*.txt ====
abundance_files = sorted(glob(os.path.join(abundance_dir, "adjusted_abundance_*.txt")))
merged_df = pd.DataFrame()

for file in abundance_files:
    sample = os.path.basename(file).replace("adjusted_abundance_", "").replace(".txt", "")
    df = pd.read_csv(file, sep="\t")

    # 保留显著项（filtered == no）
    df = df[df["filtered"] == "no"][["taxa.name", "adjusted_abundance"]]
    df.rename(columns={"adjusted_abundance": sample}, inplace=True)

    if merged_df.empty:
        merged_df = df
    else:
        merged_df = pd.merge(merged_df, df, on="taxa.name", how="outer")

# 缺失值填 0
merged_df.fillna(0, inplace=True)

# ==== 合并分类信息 ====
final_df = pd.merge(taxonomy_df, merged_df, on="taxa.name", how="right")

# ==== 保存结果 ====
final_df.to_csv(output_file, sep="\t", index=False)

print(f"✅ 合并完成，输出文件已保存到: {output_file}")
