import pandas as pd
import os

# 配置路径
base_dir = "/path/to/data2/CRC/CCDC1_TCG/FengQ_2015/Results/co_abundance_network"
groups = ["CRC", "control"]
pval_cutoff = 0.001

# 读取映射表（包含 ID、GTDBTK、fa_filename）
mapping_path = "/path/to/data2/CRC/CCDC1_TCG/FengQ_2015/Results/gtdbtk_id_mapping_with_unclassified.tsv"
mapping_df = pd.read_csv(mapping_path, sep="\t")

# 构建映射字典：fa_filename -> (ID, GTDBTK)
fa_to_id_gtdb = {}
for _, row in mapping_df.iterrows():
    fa = row["fa_filename"]
    id_ = row["ID"]
    gtdb = row["GTDBTK"]
    fa_to_id_gtdb[fa] = (id_, gtdb)

for group in groups:
    input_dir = os.path.join(base_dir, group)
    print(f"\n📂 正在处理: {group}")

    corr_path = os.path.join(input_dir, "correlation.tsv")
    pval_path = os.path.join(input_dir, "pvalues.tsv")

    corr = pd.read_csv(corr_path, sep="\t", index_col=0)
    pval = pd.read_csv(pval_path, sep="\t", index_col=0)

    edges = []
    nodes_info = {}

    for i in corr.index:
        for j in corr.columns:
            if i < j:
                p = pval.loc[i, j]
                if p <= pval_cutoff:
                    c = corr.loc[i, j]

                    id1, gtdb1 = fa_to_id_gtdb.get(i, (i, "Unclassified"))
                    id2, gtdb2 = fa_to_id_gtdb.get(j, (j, "Unclassified"))
                    node1 = id1
                    node2 = id2

                    edges.append([node1, node2, c, p])
                    nodes_info[node1] = {"ID": id1, "GTDBTK": gtdb1}
                    nodes_info[node2] = {"ID": id2, "GTDBTK": gtdb2}

    edges_df = pd.DataFrame(edges, columns=["Source", "Target", "Correlation", "Pvalue"])
    n_nodes = len(nodes_info)
    n_edges = len(edges_df)

    # 输出边文件（使用 ID 作为节点名，确保唯一）
    edge_filename = f"{group}_network_edges_({n_nodes}nodes_{n_edges}edges).tsv"
    edge_path = os.path.join(input_dir, edge_filename)
    edges_df.to_csv(edge_path, sep="\t", index=False)

    # 输出注释文件：Node（ID）、GTDBTK 分开列
    node_df = pd.DataFrame.from_dict(nodes_info, orient="index").reset_index()
    node_df.rename(columns={"index": "Node"}, inplace=True)

    anno_path = os.path.join(input_dir, f"{group}_node_annotations.tsv")
    node_df.to_csv(anno_path, sep="\t", index=False)

    print(f"✅ {group}: {n_nodes} 个节点，{n_edges} 条显著边，已输出边文件和注释文件。")
