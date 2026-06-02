import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm
from matplotlib.ticker import MaxNLocator

# === 1. 读取数据 ===
df = pd.read_csv("/path/to/data2/CRC/CCDC1_TCG/FengQ_2015/Results/co_abundance_network/2_letter_pattern_counts.tsv", sep="\t")
df = df.sort_values(by="Count", ascending=False)
df["Percent"] = df["Count"] / df["Count"].sum() * 100

labels = df["Status"]
counts = df["Count"]
y_pos = range(len(df))

# === 读取 stable_edge_node_counts.csv 文件 ===
edge_node_df = pd.read_csv("/path/to/data2/CRC/CCDC1_TCG/FengQ_2015/Results/co_abundance_network/stable_edge_node_counts.csv")
# 假设文件内容类似于：稳定边数, 节点数
stable_edge = edge_node_df.iloc[0, 0]  # 第一行第一列是稳定边数
node_count = edge_node_df.iloc[0, 1]  # 第一行第二列是节点数

# === 2. 使用颜色映射（tab20最多支持20个不同色） ===
colors = [cm.tab20(i % 20) for i in range(len(df))]

# === 3. 创建三段断轴图 ===
fig, (ax1, ax2, ax3) = plt.subplots(1, 3, sharey=True, figsize=(10, 6),
                                    gridspec_kw={'width_ratios': [1, 1, 1]},
                                    constrained_layout=True)  # 更智能布局

# === 设置每个子图的 x 范围 ===
ax1.set_xlim(0, 25) 
ax1.xaxis.set_major_locator(MaxNLocator(integer=True))
ax2.set_xlim(30, 100)
ax3.set_xlim(500, 4000)

# === 4. 绘图：每个子图只画自己范围内可见的部分，避免断轴边界出现伪影 ===
def plot_visible_segment(ax, xmin, xmax):
    visible_width = (counts.clip(lower=xmin, upper=xmax) - xmin).clip(lower=0)
    visible_mask = visible_width > 0

    if visible_mask.any():
        visible_y = [y for y, keep in zip(y_pos, visible_mask) if keep]
        visible_colors = [color for color, keep in zip(colors, visible_mask) if keep]
        ax.barh(
            visible_y,
            visible_width[visible_mask],
            left=xmin,
            color=visible_colors,
            height=0.8,
            linewidth=0,
            edgecolor="none",
        )

    ax.set_ylim(-0.5, len(df) - 0.5)

plot_visible_segment(ax1, 0, 25)
plot_visible_segment(ax2, 30, 100)
plot_visible_segment(ax3, 500, 4000)

# === 5. 添加百分比标签（支持 %, ‰, ‱）=== 
for i, (x, p) in enumerate(zip(counts, df["Percent"])): 
    if x == 0:
        continue  # 跳过 0，不画标签

    if p >= 0.1:
        label = f"{p:.2f}%"
    elif p >= 0.01:
        label = f"{p*10:.2f}‰"
    else:
        label = f"{p*100:.2f}‱"

    if x <= 25:
        ax1.text(x + 0.5, i, label, va="center", fontsize=8)
    elif x <= 100:
        ax2.text(x + 20, i, label, va="center", fontsize=8)
    else:
        ax3.text(x + 300, i, label, va="center", fontsize=8)

# === 6. Y轴标签仅左图保留 ===
ax1.set(yticks=y_pos, yticklabels=labels)
ax1.tick_params(labelsize=10)
ax1.yaxis.tick_left()
ax1.set_ylabel("Types", fontsize=18)
ax2.tick_params(axis='y', left=False)
ax3.tick_params(axis='y', left=False)

# === 7. 移除子图之间的分隔线 ===
for ax in (ax1, ax2, ax3):
    ax.spines['right'].set_visible(False)

ax2.spines['left'].set_visible(False)
ax3.spines['left'].set_visible(False)

# === 8. 添加断轴符号（斜线）=== 
d = .015
kwargs = dict(transform=ax1.transAxes, color='k', clip_on=False)
ax1.plot((1 - d, 1 + d), (-d, +d), **kwargs)
ax1.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)
kwargs.update(transform=ax2.transAxes)
ax2.plot((-d, +d), (-d, +d), **kwargs)
ax2.plot((-d, +d), (1 - d, 1 + d), **kwargs)

kwargs = dict(transform=ax2.transAxes, color='k', clip_on=False)
ax2.plot((1 - d, 1 + d), (-d, +d), **kwargs)
ax2.plot((1 - d, 1 + d), (1 - d, 1 + d), **kwargs)
kwargs.update(transform=ax3.transAxes)
ax3.plot((-d, +d), (-d, +d), **kwargs)
ax3.plot((-d, +d), (1 - d, 1 + d), **kwargs)

# === 9. 标题和标签 ===
fig.suptitle("Correlation Across Networks", fontsize=16)
ax2.set_xlabel("Number of genome pairs")
ax1.invert_yaxis()

# === 10. 在图中添加文本 ===
text = f"{stable_edge} stable genome pairs from {node_count} genomes"
fig.text(0.75, 0.2, text, ha='center', va='center', fontsize=12, color='black', fontweight='bold')

# === 11. 保存图 ===
plt.savefig("/path/to/data2/CRC/CCDC1_TCG/FengQ_2015/Results/co_abundance_network/correlation_pattern_barplot.pdf", dpi=300)
plt.show()
