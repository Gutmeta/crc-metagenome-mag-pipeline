library(dplyr)
library(stringr)
library(dendextend)

setwd("/path/to/data2/CRC/CCDC1_TCG/FengQ_2015/Results")

# === 读取边文件 ===
edges <- read.csv("C1_edges.csv", header = TRUE)
edges <- edges %>%
  mutate(Source = str_extract(name, "^[^ ]+"),
         Target = str_extract(name, "(?<=\\) )[A-Za-z0-9_.-]+")) %>%
  na.omit() %>%
  select(Source, Target, Correlation)

# === 读取注释 ===
anno <- read.table("gtdbtk_id_mapping_with_unclassified.tsv", sep = "\t", header = TRUE) # nolint
anno <- anno %>% mutate(NewName = paste(GTDBTK, ID, sep = " "))
id_to_name <- setNames(anno$NewName, anno$ID)

# === 替换名称 ===
edges$Source <- ifelse(edges$Source %in% names(id_to_name), id_to_name[edges$Source], edges$Source) # nolint: line_length_linter.
edges$Target <- ifelse(edges$Target %in% names(id_to_name), id_to_name[edges$Target], edges$Target)

# === 构建邻接矩阵 ===
nodes <- unique(c(edges$Source, edges$Target))
adj_mat <- matrix(0, length(nodes), length(nodes), dimnames = list(nodes, nodes))
for (i in 1:nrow(edges)) {
  a <- edges$Source[i]
  b <- edges$Target[i]
  c <- edges$Correlation[i]
  adj_mat[a, b] <- c
  adj_mat[b, a] <- c
}

# === 聚类 ===
dist_mat <- as.dist(1 - adj_mat)
hc <- hclust(dist_mat, method = "average")
dend <- as.dendrogram(hc)

# === 分模块 ===
clusters <- cutree(hc, k = 2)
labels_colors(dend) <- ifelse(labels(dend) %in% names(clusters[clusters == 1]), "chartreuse4", "purple3")
group_labels <- ifelse(clusters[labels(dend)] == 1, "cluster1", "cluster2")
group_colors <- ifelse(group_labels == "cluster1", "chartreuse4", "purple3")

# === 画图 ===
# === 画图（垂直颜色条 + 横向树） ===
pdf("Cluster_C1_Colored_Dendrogram_Vertical.pdf", width = 10, height = 40)
par(mar = c(2, 10, 2, 22))
plot(dend, horiz = TRUE, main = "Cluster C1 Substructure", cex = 0.5)


# 添加组标签（可选，标出cluster1/cluster2位置）
text(x = 1.1, y = mean(which(group_labels == "cluster1")), labels = "cluster1", col = "chartreuse4", xpd = TRUE, cex = 1.5)
text(x = 1.1, y = mean(which(group_labels == "cluster2")), labels = "cluster2", col = "purple3", xpd = TRUE, cex = 1.5)

dev.off()
cat("✅ 聚类图已保存为 Cluster_C1_Colored_Dendrogram_Vertical.pdf\n")

# === 生成分组表格（节点 -> cluster1/cluster2） ===
group_df <- data.frame(
  Node = labels(dend),
  Group = group_labels,
  stringsAsFactors = FALSE
)

# 保存为 CSV 文件
write.csv(group_df, "C1_clusters_table.csv", row.names = FALSE)
cat("✅ 分组表格已保存为 C1_clusters_table.csv\n")

# === WGCNA分析（识别模块）===
#library(WGCNA)

# WGCNA expects samples as rows and features as columns; transpose the adjacency matrix to form a pseudo-expression matrix.
#fake_expr <- t(adj_mat)

# 关闭WGCNA的交互式提示
#options(stringsAsFactors = FALSE)
#allowWGCNAThreads()  # Optional multithreading.

# 构建模块
#net <- blockwiseModules(
#  fake_expr,
#  power = 1,                      # Binary signed adjacency values use power 1.
#  TOMType = "unsigned",          # 不考虑正负号方向
#  minModuleSize = 10,            # Minimum module size.
# reassignThreshold = 0,
# mergeCutHeight = 0.25,
#  numericLabels = TRUE,
#  pamRespectsDendro = FALSE,
#  verbose = 3
#)

# === 提取模块颜色和标签 ===
#moduleLabels <- net$colors
#moduleColors <- labels2colors(moduleLabels)

# 添加模块信息到节点表
#module_df <- data.frame(Node = rownames(fake_expr), Module = moduleColors)

# 导出模块信息
#write.table(module_df, "Cluster_C1_WGCNA_Modules.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
#cat("✅ WGCNA模块分组已保存为 Cluster_C1_WGCNA_Modules.tsv\n")

# 可视化模块树状图
#pdf("Cluster_C1_WGCNA_Dendrogram.pdf", width = 12, height = 8)
#plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
#                    "Module colors", dendroLabels = FALSE, hang = 0.03,
#                    addGuide = TRUE, guideHang = 0.05)
#dev.off()
#cat("✅ WGCNA聚类图已保存为 Cluster_C1_WGCNA_Dendrogram.pdf\n")
