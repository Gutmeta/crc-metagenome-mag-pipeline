import glob
import pandas as pd
import csv
import os
from collections import defaultdict
from itertools import product

# === 参数设定 ===
base_dir = "/path/to/data2/CRC/CCDC1_TCG/FengQ_2015/Results/co_abundance_network"
# 去掉 adenoma，只保留两组
groups = ["CRC", "control"]
mapping_path = "/path/to/data2/CRC/CCDC1_TCG/FengQ_2015/Results/gtdbtk_id_mapping_with_unclassified.tsv"
output_edge_file = os.path.join(base_dir, "stable_network_edges.tsv")
output_node_file = os.path.join(base_dir, "stable_network_nodes.tsv")

# === 读取 GTDBTK 注释映射表：ID -> GTDBTK taxonomy
mapping_df = pd.read_csv(mapping_path, sep="\t")
id_to_gtdb = dict(zip(mapping_df["ID"], mapping_df["GTDBTK"]))

# === 累积边状态向量 v_ij，长度=组数
edge_status = defaultdict(lambda: [0] * len(groups))  # {(node1, node2): [v1, v2, ...]}

for g_idx, group in enumerate(groups):
    pattern = os.path.join(base_dir, group, f"{group}_network_edges_(*nodes_*edges).tsv")
    matched_files = glob.glob(pattern)
    if not matched_files:
        print(f"⚠️ 未找到匹配文件: {pattern}")
        continue

    edge_file = matched_files[0]
    df = pd.read_csv(edge_file, sep="\t")

    # 更稳健的 Source/Target 列识别
    lower_map = {c.lower(): c for c in df.columns}
    if "source" in lower_map and "target" in lower_map:
        source_col, target_col = lower_map["source"], lower_map["target"]
    else:
        raise ValueError(f"❌ 边表缺失 Source/Target 列: {df.columns.tolist()}")

    if "correlation" in lower_map:
        corr_col = lower_map["correlation"]
    else:
        raise ValueError(f"❌ 边表缺失 Correlation 列: {df.columns.tolist()}")

    for _, row in df.iterrows():
        i, j = sorted([row[source_col], row[target_col]])
        corr = row[corr_col]
        if corr > 0:
            edge_status[(i, j)][g_idx] = 1
        elif corr < 0:
            edge_status[(i, j)][g_idx] = -1
        # corr == 0 或缺失时保持 0（U）

# === 提取稳定边（在所有组中都一致为 +1 或 -1）
stable_edges = []
node_set = set()

for (i, j), vec in edge_status.items():
    total = sum(vec)
    if total == len(groups):      # 全正相关（两组 -> +2）
        stable_edges.append([i, j, 1])
        node_set.update([i, j])
    elif total == -len(groups):   # 全负相关（两组 -> -2）
        stable_edges.append([i, j, -1])
        node_set.update([i, j])

print(f"📊 共提取到稳定边数: {len(stable_edges)}，涉及节点数: {len(node_set)}")
filename = 'stable_edge_node_counts.csv'
file_path = os.path.join(base_dir, filename)
with open(file_path, 'w', newline='', encoding='utf-8') as csvfile:
    writer = csv.writer(csvfile)
    
    # 写入表头
    writer.writerow(["稳定边数", "节点数"])
    
    # 写入数据行
    writer.writerow([len(stable_edges), len(node_set)])

# === 统计字母组合的个数（按组数自适应；两组即两字母）
status_count = defaultdict(int)

def vec_to_label(vec):
    # 1 -> P，-1 -> N，0 -> U
    return ''.join('P' if v == 1 else 'N' if v == -1 else 'U' for v in vec)

# 统计实际出现的组合
for vec in edge_status.values():
    status_count[vec_to_label(vec)] += 1

# 构造所有可能的组合（两组 -> PP/PN/NP/NN/PU/UP/NU/UN/UU）
all_combinations = [''.join(p) for p in product("PUN", repeat=len(groups))]

# 先统计“非全 U”组合的总数
non_all_u_total = sum(v for k, v in status_count.items() if set(k) != {"U"})

# 用 mapping_df 中的 ID 数量估算“所有可能边数”
total_nodes = mapping_df["ID"].nunique()
total_possible_edges = total_nodes * (total_nodes - 1) // 2

# 组装完整计数，并把“全 U（例如两组时为 'UU'）”修正为：总可能边数 − 非全 U
full_status_counts = []
all_u_label = "U" * len(groups)
for comb in all_combinations:
    if comb == all_u_label:
        count_val = total_possible_edges - non_all_u_total
    else:
        count_val = status_count.get(comb, 0)
    full_status_counts.append({"Status": comb, "Count": int(count_val)})

# 输出组合统计
status_df = pd.DataFrame(full_status_counts).sort_values(by="Count", ascending=False)
output_path = os.path.join(base_dir, f"{len(groups)}_letter_pattern_counts.tsv")
status_df.to_csv(output_path, sep="\t", index=False)
print(f"✅ 组合统计表已保存为: {output_path}")

# === 构建稳定边表
edge_df = pd.DataFrame(stable_edges, columns=["Source", "Target", "Correlation"])
edge_df.to_csv(output_edge_file, sep="\t", index=False)
print(f"✅ 边文件已保存: {output_edge_file}")

# === 构建节点注释表（ID -> GTDBTK）
node_records = []
for node_id in node_set:
    gtdb = id_to_gtdb.get(node_id, "Unclassified")
    node_records.append([node_id, gtdb])

node_df = pd.DataFrame(node_records, columns=["Node", "GTDBTK"]).drop_duplicates()
node_df.to_csv(output_node_file, sep="\t", index=False)
print(f"✅ 节点注释文件已保存: {output_node_file}")
