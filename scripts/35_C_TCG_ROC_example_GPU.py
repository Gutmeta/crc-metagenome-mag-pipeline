import os
import pandas as pd
import numpy as np
import cupy as cp
from glob import glob

from sklearn.model_selection import LeaveOneOut
from sklearn.metrics import roc_curve, auc, accuracy_score, precision_score, recall_score, f1_score
from sklearn.calibration import calibration_curve
from sklearn.metrics import confusion_matrix

import matplotlib.pyplot as plt
from xgboost import XGBClassifier


# =========================
# 0) 统一路径配置
# =========================
DATASET_NAME = "FengQ_2015"

ROOT_DATA5 = "/path/to/data2/CRC/CCDC2"
ROOT_DATA2 = "/path/to/data1"

DATA_DIR = os.path.join(ROOT_DATA5, DATASET_NAME, "DiTASiC", "abundance")
SEQKIT_STATS_FILE = os.path.join(ROOT_DATA2, DATASET_NAME, "seqkit_stats.txt")

# Group metadata file, resolved relative to the working directory.
GROUP_FILE = f"{DATASET_NAME}_CRC_Group.txt"

OUT_DIR = os.path.join(ROOT_DATA5, "ML_results")
os.makedirs(OUT_DIR, exist_ok=True)

# 输出文件前缀（统一管理）
OUT_PREFIX = os.path.join(OUT_DIR, DATASET_NAME)


# =========================
# 1) 读取分组信息
# =========================
group_df = pd.read_csv(GROUP_FILE, sep="\t", header=0)
group_df.columns = ['Sample', 'Group']
group_df['label'] = group_df['Group'].fillna('').astype(str).apply(lambda x: 1 if x.strip() == 'CRC' else 0)
group_df = group_df.set_index('Sample')


# =========================
# 2) 读取 seqkit_stats.txt
# =========================
seqkit_df = pd.read_csv(SEQKIT_STATS_FILE, sep=r'\s+', header=0)

# 提取样本ID与reads数
seqkit_df['SampleID'] = seqkit_df['file'].str.split('/').str[0]
seqkit_df['num_seqs'] = seqkit_df['num_seqs'].astype(str).str.replace(',', '').astype(int)
seqkit_df['num_reads'] = seqkit_df['num_seqs'].astype(int)

total_reads_per_sample = seqkit_df.groupby('SampleID')['num_reads'].sum()


# =========================
# 3) 读取 abundance 文件并合并成矩阵（剔除空样本/缺reads样本）
# =========================
file_list = glob(os.path.join(DATA_DIR, "abundance_*.txt"))

# 诊断：group 里有哪些样本根本没有 abundance 文件
file_samples = set(
    os.path.basename(p).replace("abundance_", "").replace(".txt", "")
    for p in file_list
)
group_samples = set(group_df.index)
missing_abundance_files = sorted(group_samples - file_samples)
print(f"[INFO] group里但没有abundance文件的样本数: {len(missing_abundance_files)}")
print(f"[INFO] 示例(前20): {missing_abundance_files[:20]}")

sample_data = []
kept_samples = []

skipped_not_in_group = []
skipped_no_total_reads = []
skipped_empty_after_pval = []

for file_path in file_list:
    sample_id = os.path.basename(file_path).replace("abundance_", "").replace(".txt", "")

    # 1) 不在分组表里的文件直接跳过
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

    # 4) Drop samples with no taxa after p-value filtering to avoid all-zero rows.
    if df.empty:
        skipped_empty_after_pval.append(sample_id)
        continue

    s = df.set_index("taxa.name")["count.estimate"]
    s.name = sample_id

    sample_data.append(s)
    kept_samples.append(sample_id)

print(f"[INFO] kept (used) samples: {len(kept_samples)}")
print(f"[INFO] skipped_not_in_group: {len(skipped_not_in_group)}")
print(f"[INFO] skipped_no_total_reads: {len(skipped_no_total_reads)}")
print(f"[INFO] skipped_empty_after_pval: {len(skipped_empty_after_pval)}")

# 合并矩阵
abundance_matrix = pd.DataFrame(sample_data).fillna(0)

# reads 严格按 abundance_matrix 的样本顺序对齐
reads_aligned = total_reads_per_sample.loc[abundance_matrix.index].astype(float)

# 额外防御：reads<=0 的样本剔除
valid_reads_mask = reads_aligned > 0
if (~valid_reads_mask).any():
    bad_samples = reads_aligned.index[~valid_reads_mask].tolist()
    print(f"[WARN] reads<=0 的样本将被剔除: {len(bad_samples)} (前20: {bad_samples[:20]})")
    abundance_matrix = abundance_matrix.loc[valid_reads_mask]
    reads_aligned = reads_aligned.loc[valid_reads_mask]

# Normalize to relative abundance without masking invalid values.
abundance_rel = abundance_matrix.div(reads_aligned, axis=0)

# 防御性过滤：剔除全 NaN / 全 0 行
abundance_rel = abundance_rel.dropna(axis=0, how="all")
abundance_rel = abundance_rel.loc[(abundance_rel.fillna(0) != 0).any(axis=1)]

# Remove all-zero columns to reduce unused features.
abundance_rel = abundance_rel.loc[:, (abundance_rel.fillna(0) != 0).any(axis=0)]

# 对齐标签（最终进入模型样本）
abundance_rel = abundance_rel.loc[group_df.index.intersection(abundance_rel.index)]
labels = group_df.loc[abundance_rel.index, 'label']

# 最终 X / y
X = abundance_rel.to_numpy(dtype=np.float32)
y = labels.to_numpy(dtype=int)

# 特征名
feature_names = abundance_rel.columns.tolist()
n_features = len(feature_names)

print(f"[INFO] Dataset: {DATASET_NAME}")
print(f"[INFO] Samples: {X.shape[0]}, Features: {X.shape[1]}")



# =========================
# 4) Leave-One-Out 交叉验证 + 记录特征重要性
# =========================
loo = LeaveOneOut()

y_true, y_scores = [], []
train_accuracies, val_accuracies = [], []

# Accumulate mean feature importance across LOOCV folds.
importance_sum = np.zeros(n_features, dtype=float)
importance_sq_sum = np.zeros(n_features, dtype=float)
n_folds = 0

for train_idx, test_idx in loo.split(X, y):
    X_train, X_test = X[train_idx], X[test_idx]
    y_train, y_test = y[train_idx], y[test_idx]

    X_train_gpu = cp.asarray(X_train, dtype=cp.float32)
    y_train_gpu = cp.asarray(y_train)
    X_test_gpu  = cp.asarray(X_test,  dtype=cp.float32)

    clf = XGBClassifier(
        eval_metric='logloss',
        random_state=42,
        tree_method='hist',
        device="cuda:0",
        # Fix the importance type to gain for consistency across folds.
        importance_type="gain"
    )

    # 训练
    clf.fit(X_train_gpu, y_train_gpu)

    # Collect feature importance for the current fold.
    fold_imp = np.asarray(clf.feature_importances_, dtype=float)
    if fold_imp.shape[0] != n_features:
        raise ValueError(
            f"Feature importance length mismatch: got {fold_imp.shape[0]}, expected {n_features}"
        )

    importance_sum += fold_imp
    importance_sq_sum += fold_imp ** 2
    n_folds += 1

    # 训练集/验证集准确率
    train_pred = cp.asnumpy(clf.predict(X_train_gpu))
    test_pred  = cp.asnumpy(clf.predict(X_test_gpu))

    train_accuracy = accuracy_score(y_train, train_pred)
    val_accuracy   = accuracy_score(y_test, test_pred)

    train_accuracies.append(train_accuracy)
    val_accuracies.append(val_accuracy)

    # 预测概率
    prob = cp.asnumpy(clf.predict_proba(X_test_gpu)[:, 1])[0]
    y_scores.append(float(prob))
    y_true.append(int(y_test[0]))


# =========================
# 5) 输出：LOOCV 平均特征重要性（均值 + 标准差 + 排名）
# =========================
importance_mean = importance_sum / max(n_folds, 1)

# 方差 = E[x^2] - (E[x])^2
importance_var = (importance_sq_sum / max(n_folds, 1)) - (importance_mean ** 2)
importance_var = np.maximum(importance_var, 0)
importance_std = np.sqrt(importance_var)

importance_df = pd.DataFrame({
    "Feature": feature_names,
    "MeanImportance": importance_mean,
    "StdImportance": importance_std
})

# 排名：均值越大越重要（rank=1 最重要）
importance_df["ImportanceRank"] = importance_df["MeanImportance"].rank(ascending=False, method="min").astype(int)
importance_df = importance_df.sort_values("MeanImportance", ascending=False)

importance_outfile = f"{OUT_PREFIX}_feature_mean_importance.csv"
importance_df.to_csv(importance_outfile, index=False)

print(f"[INFO] Saved mean feature importance to: {importance_outfile}")
print("[INFO] Top 20 features by mean importance:")
print(importance_df.head(20))


# =========================
# 6) ROC / AUC + 95% CI
# =========================
fpr, tpr, _ = roc_curve(y_true, y_scores)
roc_auc = auc(fpr, tpr)

# Bootstrapping 计算 95% CI
n_bootstraps = 1000
rng = np.random.RandomState(42)
bootstrapped_scores = []

y_true_arr = np.array(y_true)
y_scores_arr = np.array(y_scores)

for i in range(n_bootstraps):
    indices = rng.randint(0, len(y_scores_arr), len(y_scores_arr))
    if len(np.unique(y_true_arr[indices])) < 2:
        continue
    fpr_b, tpr_b, _ = roc_curve(y_true_arr[indices], y_scores_arr[indices])
    score = auc(fpr_b, tpr_b)
    bootstrapped_scores.append(score)

ci_lower = np.percentile(bootstrapped_scores, 2.5)
ci_upper = np.percentile(bootstrapped_scores, 97.5)

# 绘制 ROC 曲线
plt.figure(figsize=(6, 6))
plt.plot(fpr, tpr, label=f"AUC = {roc_auc:.2f}\n95% CI: {ci_lower:.2f}-{ci_upper:.2f}", linewidth=2)
plt.plot([0, 1], [0, 1], 'k--', lw=1)

plt.xlabel('1 - Specificity')
plt.ylabel('Sensitivity')
plt.title(f'CRC {DATASET_NAME}')
plt.legend(loc='lower right')
plt.grid(False)
plt.tight_layout()

roc_pdf = f"{OUT_PREFIX}_roc_curve.pdf"
plt.savefig(roc_pdf)
plt.close()


# =========================
# 7) 校准曲线
# =========================
prob_true, prob_pred = calibration_curve(y_true, y_scores, n_bins=10)
plt.figure(figsize=(6, 6))
plt.plot(prob_pred, prob_true, marker='o', label="Calibration curve")
plt.plot([0, 1], [0, 1], linestyle='--', label="Perfectly calibrated")
plt.xlabel("Mean predicted probability")
plt.ylabel("Fraction of positives")
plt.title("Calibration Curve")
plt.legend(loc="lower right")
plt.tight_layout()

cal_pdf = f"{OUT_PREFIX}_calibration_curve.pdf"
plt.savefig(cal_pdf)
plt.close()


# =========================
# 8) DCA 净效益
# =========================
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

thresh_group = np.arange(0, 1, 0.01)
net_benefit_model = calculate_net_benefit_model(thresh_group, np.array(y_scores), np.array(y_true))
net_benefit_all = calculate_net_benefit_all(thresh_group, np.array(y_true))

fig, ax = plt.subplots(figsize=(6, 6))
ax.plot(thresh_group, net_benefit_model, color='crimson', label='Model')
ax.plot(thresh_group, net_benefit_all, color='black', label='Treat all')
ax.plot((0, 1), (0, 0), color='black', linestyle=':', label='Treat none')

y2 = np.maximum(net_benefit_all, 0)
y1 = np.maximum(net_benefit_model, y2)
ax.fill_between(thresh_group, y1, y2, color='crimson', alpha=0.2)

ax.set_xlim(0, 1)
ax.set_ylim(min(net_benefit_model.max() - 0.15, -0.18), net_benefit_model.max() + 0.15)
ax.set_xlabel('Threshold Probability', fontdict={'family': 'Times New Roman', 'fontsize': 15})
ax.set_ylabel('Net Benefit', fontdict={'family': 'Times New Roman', 'fontsize': 15})
ax.grid('major')
ax.spines['right'].set_color((0.8, 0.8, 0.8))
ax.spines['top'].set_color((0.8, 0.8, 0.8))
ax.legend(loc='upper right')
plt.tight_layout()

dca_pdf = f"{OUT_PREFIX}_dca_curve.pdf"
plt.savefig(dca_pdf)
plt.close()


# =========================
# 9) 输出每个样本预测结果 + 指标
# =========================
sample_names = abundance_rel.index.tolist()

results_df = pd.DataFrame({
    "Sample": sample_names,
    "True Label": y_true,
    "Predicted Probability": y_scores,
})

results_df["Predicted Label"] = (results_df["Predicted Probability"] >= 0.5).astype(int)
results_df["Correct"] = results_df["True Label"] == results_df["Predicted Label"]

accuracy = accuracy_score(y_true, (np.array(y_scores) >= 0.5).astype(int))
precision = precision_score(y_true, (np.array(y_scores) >= 0.5).astype(int))
recall = recall_score(y_true, (np.array(y_scores) >= 0.5).astype(int))
f1 = f1_score(y_true, (np.array(y_scores) >= 0.5).astype(int))

metrics_df = pd.DataFrame({
    "Metric": ["Accuracy", "Precision", "Recall", "F1 Score"],
    "Score": [accuracy, precision, recall, f1]
})

metrics_csv = f"{OUT_PREFIX}_metrics.csv"
pred_csv = f"{OUT_PREFIX}_prediction_results.csv"

metrics_df.to_csv(metrics_csv, index=False)
results_df.to_csv(pred_csv, index=False)

print(metrics_df)

print("❌ 预测错误的样本：")
print(results_df[~results_df["Correct"]])

print(f"[INFO] Saved ROC: {roc_pdf}")
print(f"[INFO] Saved Calibration: {cal_pdf}")
print(f"[INFO] Saved DCA: {dca_pdf}")
print(f"[INFO] Saved metrics: {metrics_csv}")
print(f"[INFO] Saved predictions: {pred_csv}")
