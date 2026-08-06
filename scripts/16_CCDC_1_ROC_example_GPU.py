import os
from pathlib import Path
import pandas as pd
import numpy as np
import cupy as cp
from glob import glob
from sklearn.model_selection import LeaveOneOut
from sklearn.metrics import roc_curve, auc, accuracy_score, precision_score, recall_score, f1_score
from sklearn.calibration import calibration_curve
import matplotlib.pyplot as plt
from matplotlib import font_manager
from xgboost import XGBClassifier, DMatrix
from sklearn.metrics import confusion_matrix

def configure_plot_font():
    font_candidates = [
        Path.home() / ".local/share/fonts/arial-corefonts/Arial.TTF",
        Path("/usr/share/fonts/truetype/msttcorefonts/Arial.ttf"),
        Path("/usr/share/fonts/truetype/msttcorefonts/Arial.TTF"),
        Path("/usr/share/fonts/truetype/msttcorefonts/arial.ttf"),
    ]

    font_family = "DejaVu Sans"
    for font_path in font_candidates:
        if font_path.exists():
            font_manager.fontManager.addfont(str(font_path))
            font_family = font_manager.FontProperties(fname=str(font_path)).get_name()
            break

    plt.rcParams["font.family"] = font_family
    plt.rcParams["font.sans-serif"] = [font_family]


configure_plot_font()

# === 设置路径 ===
data_dir = "/path/to/data2/CRC/CCDC1/FengQ_2015/DiTASiC/abundance/"
group_file = "FengQ_2015_CRC_Group.txt"
os.makedirs("/path/to/data2/CRC/CCDC1/ML_results/", exist_ok=True)

# === 读取分组信息 ===
group_df = pd.read_csv(group_file, sep="\t", header=0)
group_df.columns = ['Sample', 'Group']
group_df['label'] = group_df['Group'].fillna('').astype(str).apply(lambda x: 1 if x.strip() == 'CRC' else 0)
group_df = group_df.set_index('Sample')

# === 读取 seqkit_stats.txt 文件 ===
seqkit_df = pd.read_csv("/path/to/data1/FengQ_2015/seqkit_stats.txt", sep=r'\s+', header=0)

# 提取样本ID和总reads数
seqkit_df['SampleID'] = seqkit_df['file'].str.split('/').str[0]  # 从文件路径中提取样本ID
seqkit_df['num_seqs'] = seqkit_df['num_seqs'].str.replace(',', '').astype(int)
seqkit_df['num_reads'] = seqkit_df['num_seqs'].astype(int)

# === 读取 abundance 文件 ===
file_list = glob(os.path.join(data_dir, "abundance_*.txt"))
sample_data = []
kept_samples = []

skipped_not_in_group = []
skipped_empty_after_pval = []
skipped_no_total_reads = []

# 先准备 total_reads_per_sample
total_reads_per_sample = seqkit_df.groupby('SampleID')['num_reads'].sum()

for file_path in file_list:
    sample_id = os.path.basename(file_path).replace("abundance_", "").replace(".txt", "")

    # 1) 不在分组表里的样本直接跳过
    if sample_id not in group_df.index:
        skipped_not_in_group.append(sample_id)
        continue

    # 2) 没有对应 total reads 的样本跳过（否则会除出 NaN）
    if sample_id not in total_reads_per_sample.index:
        skipped_no_total_reads.append(sample_id)
        continue

    df = pd.read_csv(file_path, sep="\t")

    # 3) 过滤显著 taxa
    df = df[df["raw.pval"] <= 0.05]

    # 4) Skip samples with no retained taxa to avoid all-zero rows.
    if df.empty:
        skipped_empty_after_pval.append(sample_id)
        continue

    s = df.set_index("taxa.name")["count.estimate"]
    s.name = sample_id

    sample_data.append(s)
    kept_samples.append(sample_id)

print("✅ kept (used) samples:", len(kept_samples))
print("❌ skipped_not_in_group:", len(skipped_not_in_group))
print("❌ skipped_no_total_reads:", len(skipped_no_total_reads))
print("❌ skipped_empty_after_pval:", len(skipped_empty_after_pval))

# === 合并丰度矩阵 ===
abundance_matrix = pd.DataFrame(sample_data).fillna(0)

# === 标准化为相对丰度（按样本总reads）===
# 这里用 abundance_matrix.index 来取对应 reads，确保顺序严格一致
reads_aligned = total_reads_per_sample.loc[abundance_matrix.index]

abundance_rel = abundance_matrix.div(reads_aligned, axis=0)

# Remove all-NaN and all-zero rows as a defensive validation step.
abundance_rel = abundance_rel.dropna(axis=0, how="all")
abundance_rel = abundance_rel.loc[(abundance_rel.fillna(0) != 0).any(axis=1)]

# === 对齐标签 ===
abundance_rel = abundance_rel.loc[group_df.index.intersection(abundance_rel.index)]
labels = group_df.loc[abundance_rel.index, 'label']

X = abundance_rel.values
y = labels.values

print("最终进入模型的样本数:", X.shape[0])
print("特征数(taxa):", X.shape[1])


# === Leave-One-Out 交叉验证 ===
loo = LeaveOneOut()
y_true, y_scores = [], []
train_accuracies, val_accuracies = [], []  # 用于记录训练和验证集的准确率

for train_idx, test_idx in loo.split(X):
    X_train, X_test = X[train_idx], X[test_idx]
    y_train, y_test = y[train_idx], y[test_idx]

    X_train_gpu = cp.asarray(X_train, dtype=cp.float32)
    y_train_gpu = cp.asarray(y_train)
    X_test_gpu  = cp.asarray(X_test,  dtype=cp.float32)

    # Enable GPU training with an explicit device identifier.
    clf = XGBClassifier(
        eval_metric='logloss', 
        random_state=42, 
        tree_method='hist',  # 使用 GPU 加速的树算法
        device="cuda:0"  # 显式选择使用第一个 GPU
    )
    
    # 训练模型
    clf.fit(X_train_gpu, y_train_gpu)
    
    # 获取训练集和验证集的准确率
    train_pred = cp.asnumpy(clf.predict(X_train_gpu))
    test_pred  = cp.asnumpy(clf.predict(X_test_gpu))

    train_accuracy = accuracy_score(y_train, train_pred)
    val_accuracy   = accuracy_score(y_test, test_pred)

    train_accuracies.append(train_accuracy)
    val_accuracies.append(val_accuracy)
    
    # 使用GPU进行预测
    prob = cp.asnumpy(clf.predict_proba(X_test_gpu)[:, 1])[0]
    y_scores.append(float(prob))
    y_true.append(int(y_test[0]))

# === 计算 ROC 和 AUC ===
fpr, tpr, _ = roc_curve(y_true, y_scores)
roc_auc = auc(fpr, tpr)

# === Bootstrapping 计算 95% CI ===
n_bootstraps = 1000
rng = np.random.RandomState(42)
bootstrapped_scores = []
bootstrapped_tprs = []
fpr_grid = np.linspace(0, 1, 201)

y_true_arr = np.array(y_true)
y_scores_arr = np.array(y_scores)
y_pred_arr = (y_scores_arr >= 0.5).astype(int)

for i in range(n_bootstraps):
    indices = rng.randint(0, len(y_scores_arr), len(y_scores_arr))
    if len(np.unique(y_true_arr[indices])) < 2:
        continue

    fpr_b, tpr_b, _ = roc_curve(y_true_arr[indices], y_scores_arr[indices])
    score = auc(fpr_b, tpr_b)
    bootstrapped_scores.append(score)

    interp_tpr = np.interp(fpr_grid, fpr_b, tpr_b)
    interp_tpr[0] = 0.0
    interp_tpr[-1] = 1.0
    bootstrapped_tprs.append(interp_tpr)

ci_lower = np.percentile(bootstrapped_scores, 2.5)
ci_upper = np.percentile(bootstrapped_scores, 97.5)

tpr_ci_lower = np.percentile(bootstrapped_tprs, 2.5, axis=0)
tpr_ci_upper = np.percentile(bootstrapped_tprs, 97.5, axis=0)

# === 计算分类指标（用于图上展示和导出） ===
accuracy = accuracy_score(y_true_arr, y_pred_arr)
precision = precision_score(y_true_arr, y_pred_arr, zero_division=0)
recall = recall_score(y_true_arr, y_pred_arr, zero_division=0)
f1 = f1_score(y_true_arr, y_pred_arr, zero_division=0)
tn, fp, fn, tp = confusion_matrix(y_true_arr, y_pred_arr).ravel()
specificity = tn / (tn + fp) if (tn + fp) > 0 else 0.0

# === 绘制 ROC 曲线 ===
roc_color = "#ff5a36"
roc_fill_color = "#ffd9cf"

fig, ax = plt.subplots(figsize=(7.4, 7.0))
ax.set_facecolor("#fffdfb")

ax.fill_between(
    fpr_grid * 100,
    tpr_ci_lower * 100,
    tpr_ci_upper * 100,
    step='post',
    color=roc_fill_color,
    alpha=0.45,
    linewidth=0
)
ax.step(fpr * 100, tpr * 100, where='post', color=roc_color, linewidth=2.8)
ax.plot([0, 100], [0, 100], linestyle='--', color='#cfcfcf', linewidth=1.2)

ax.set_xlim(-5, 105)
ax.set_ylim(-5, 105)
ax.set_xticks([0, 25, 50, 75, 100])
ax.set_yticks([0, 25, 50, 75, 100])
ax.set_xlabel('False Positive Percentage', fontsize=21)
ax.set_ylabel('True Positive Percentage', fontsize=21)
ax.set_title('CRC FengQ_2015', fontsize=20, pad=18)

for tick_label in ax.get_xticklabels() + ax.get_yticklabels():
    tick_label.set_fontsize(16)

ax.grid(False)
for spine in ax.spines.values():
    spine.set_linewidth(1.2)
    spine.set_color('#bdbdbd')

metric_text = "\n".join([
    f"AUROC: {roc_auc * 100:.2f}%",
    f"95% CI: {ci_lower * 100:.2f}-{ci_upper * 100:.2f}%",
    f"Accuracy: {accuracy * 100:.2f}%",
    f"Precision: {precision * 100:.2f}%",
    f"Recall: {recall * 100:.2f}%",
    f"Specificity: {specificity * 100:.2f}%",
    f"F1 score: {f1 * 100:.2f}%"
])

ax.text(
    0.98,
    0.07,
    metric_text,
    transform=ax.transAxes,
    ha='right',
    va='bottom',
    fontsize=13.5,
    color=roc_color,
    linespacing=1.35,
    bbox=dict(boxstyle='round,pad=0.4', facecolor='white', edgecolor='none', alpha=0.8)
)

plt.tight_layout()

# 保存 ROC 曲线为 PDF
plt.savefig("/path/to/data2/CRC/CCDC1/ML_results/FengQ_2015_roc_curve.pdf", bbox_inches='tight')
plt.close()

# === 校准曲线 ===
prob_true, prob_pred = calibration_curve(y_true, y_scores, n_bins=10)
plt.figure(figsize=(6, 6))
plt.plot(prob_pred, prob_true, marker='o', label="Calibration curve")
plt.plot([0, 1], [0, 1], linestyle='--', label="Perfectly calibrated")
plt.xlabel("Mean predicted probability")
plt.ylabel("Fraction of positives")
plt.title("Calibration Curve")
plt.legend(loc="lower right")
plt.tight_layout()

# 保存校准曲线为 PDF
plt.savefig("/path/to/data2/CRC/CCDC1/ML_results/FengQ_2015_calibration_curve.pdf")
plt.close()

# === 计算 DCA 净效益 ===
def calculate_net_benefit_model(thresh_group, y_pred_score, y_label):
    net_benefit_model = np.array([])  
    for thresh in thresh_group:  
        y_pred_label = y_pred_score > thresh  
        tn, fp, fn, tp = confusion_matrix(y_label, y_pred_label).ravel()  
        n = len(y_label)  
        net_benefit = (tp / n) - (fp / n) * (thresh / (1 - thresh))  
        net_benefit_model = np.append(net_benefit_model, net_benefit)  
    return net_benefit_model

def calculate_net_benefit_all(thresh_group, y_label):
    net_benefit_all = np.array([])  
    tn, fp, fn, tp = confusion_matrix(y_label, y_label).ravel()  
    total = tp + tn  
    for thresh in thresh_group:
        net_benefit = (tp / total) - (tn / total) * (thresh / (1 - thresh))  
        net_benefit_all = np.append(net_benefit_all, net_benefit)  
    return net_benefit_all

# 生成阈值范围
thresh_group = np.arange(0, 1, 0.01)
net_benefit_model = calculate_net_benefit_model(thresh_group, np.array(y_scores), np.array(y_true))
net_benefit_all = calculate_net_benefit_all(thresh_group, np.array(y_true))

# === 绘制 DCA 曲线 ===
fig, ax = plt.subplots(figsize=(6, 6))
# 绘制决策曲线分析（DCA）
ax.plot(thresh_group, net_benefit_model, color='crimson', label='Model')
ax.plot(thresh_group, net_benefit_all, color='black', label='Treat all')
ax.plot((0, 1), (0, 0), color='black', linestyle=':', label='Treat none')

y2 = np.maximum(net_benefit_all, 0)
y1 = np.maximum(net_benefit_model, y2)
ax.fill_between(thresh_group, y1, y2, color='crimson', alpha=0.2)

ax.set_xlim(0, 1)
ax.set_ylim(min(net_benefit_model.max() - 0.15, -0.18), net_benefit_model.max() + 0.15)
ax.set_xlabel('Threshold Probability', fontsize=15)
ax.set_ylabel('Net Benefit', fontsize=15)
ax.grid('major')
ax.spines['right'].set_color((0.8, 0.8, 0.8))
ax.spines['top'].set_color((0.8, 0.8, 0.8))
ax.legend(loc='upper right')
plt.tight_layout()

# 保存 DCA 曲线为 PDF
plt.savefig("/path/to/data2/CRC/CCDC1/ML_results/FengQ_2015_dca_curve.pdf")
plt.close()



# === 输出每个样本预测结果 ===
sample_names = abundance_rel.index.tolist()

# 构建 DataFrame
results_df = pd.DataFrame({
    "Sample": sample_names,
    "True Label": y_true,
    "Predicted Probability": y_scores,
})

# 判断预测是否正确（阈值默认 0.5）
results_df["Predicted Label"] = (results_df["Predicted Probability"] >= 0.5).astype(int)
results_df["Correct"] = results_df["True Label"] == results_df["Predicted Label"]

# 创建表格显示这些指标
metrics_df = pd.DataFrame({
    "Metric": ["Accuracy", "Precision", "Recall", "Specificity", "F1 Score", "AUROC", "AUROC 95% CI Lower", "AUROC 95% CI Upper"],
    "Score": [accuracy, precision, recall, specificity, f1, roc_auc, ci_lower, ci_upper]
})

# 保存指标表格为 CSV
metrics_df.to_csv("/path/to/data2/CRC/CCDC1/ML_results/FengQ_2015_metrics.csv", index=False)

# 输出指标表格
print(metrics_df)

# 保存预测结果为文件
results_df.to_csv("/path/to/data2/CRC/CCDC1/ML_results/FengQ_2015_prediction_results.csv", index=False)

# 打印错误的样本
print("❌ 预测错误的样本：")
print(results_df[~results_df["Correct"]])
