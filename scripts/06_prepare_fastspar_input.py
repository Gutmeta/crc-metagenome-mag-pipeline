import pandas as pd
import os

# 输入路径
abundance_file = "/path/to/storage/data3/CRC/YachidaS_2019/Results/merged_abundance_with_taxonomy.tsv"
group_file = "/path/to/crc-metagenome-mag-pipeline/Yachidas_2019_CRC_Group.txt"
output_base = "/path/to/storage/data3/CRC/YachidaS_2019/Results/co_abundance_network"

os.makedirs(output_base, exist_ok=True)

# 读取数据
abundance_df = pd.read_csv(abundance_file, sep="\t")
group_df = pd.read_csv(group_file, sep="\t")

# 样本名映射
sample_to_group = dict(zip(group_df["Sample"], group_df["Group"]))

# 所有样本列
sample_columns = [col for col in abundance_df.columns if col not in ["taxa.name", "taxonomy"]]

# 标准化样本名并分类
sample_mapping = {}
for col in sample_columns:
    sample_id = col.split(".")[0]
    if sample_id in sample_to_group:
        sample_mapping[col] = (sample_id, sample_to_group[sample_id])

# 分组采样列
grouped_samples = {"CRC": [], "control": []}
for col, (sid, grp) in sample_mapping.items():
    if grp in grouped_samples:
        grouped_samples[grp].append(col)

# 1️⃣ 每组筛选出现率 ≥ 0.75 的 MAG，保存列表
group_to_core_MAGs = {}

for grp, cols in grouped_samples.items():
    if not cols:
        continue
    sub_abundance = abundance_df[["taxa.name"] + cols].copy()
    presence_ratio = (sub_abundance[cols] > 0).sum(axis=1) / len(cols)
    core_MAGs = sub_abundance.loc[presence_ratio >= 0.75, "taxa.name"]
    group_to_core_MAGs[grp] = set(core_MAGs)
    print(f"✅ {grp}: 筛出核心 MAG 数量: {len(core_MAGs)}")

# 2️⃣ 取交集：核心 MAG
core_MAG_intersection = set.intersection(*group_to_core_MAGs.values())
print(f"\n🧬 所有组共有的核心 MAG 数量: {len(core_MAG_intersection)}")

# 3️⃣ 为每组生成 FastSpar 输入文件，统一使用 core_MAG_intersection
for grp, cols in grouped_samples.items():
    if not cols:
        continue

    grp_dir = os.path.join(output_base, grp)
    os.makedirs(grp_dir, exist_ok=True)

    sub_abundance = abundance_df[["taxa.name", "taxonomy"] + cols].copy()
    sub_abundance = sub_abundance[sub_abundance["taxa.name"].isin(core_MAG_intersection)].reset_index(drop=True)

    # 保存带 taxonomy 的文件
    sub_abundance.to_csv(f"{grp_dir}/filtered_with_taxonomy.tsv", sep="\t", index=False)

    # FastSpar 输入（OTU × 样本格式）
    fastspar_input = sub_abundance.drop(columns=["taxonomy"])
    fastspar_input = fastspar_input.rename(columns={"taxa.name": "#OTU ID"})
    fastspar_input.to_csv(f"{grp_dir}/fastspar_input.tsv", sep="\t", index=False)

    print(f"📁 {grp} 输出完成。样本数: {len(cols)}，MAG数: {fastspar_input.shape[0]}")

print("\n🎉 所有组处理完成，核心 MAG 统一，可进行可比网络构建。")
