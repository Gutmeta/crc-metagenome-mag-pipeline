import os
import sys
import importlib.util
from glob import glob
import numpy as np
import pandas as pd

MPLCONFIGDIR = "/tmp/codex-matplotlib"
os.makedirs(MPLCONFIGDIR, exist_ok=True)
os.environ.setdefault("MPLCONFIGDIR", MPLCONFIGDIR)

from sklearn.model_selection import LeaveOneOut
from sklearn.metrics import roc_auc_score
from scipy.stats import rankdata

CRC_PYTHON = "/path/to/conda/envs/CRC/bin/python"
if importlib.util.find_spec("xgboost") is None:
    if (
        __name__ == "__main__"
        and os.path.exists(CRC_PYTHON)
        and os.path.realpath(sys.executable) != os.path.realpath(CRC_PYTHON)
    ):
        print(f"[INFO] xgboost not found in {sys.executable}, retrying with {CRC_PYTHON}")
        os.execv(CRC_PYTHON, [CRC_PYTHON, *sys.argv])
    raise ModuleNotFoundError(
        "xgboost is not installed in the current interpreter. "
        f"Please use {CRC_PYTHON} or install xgboost in {sys.executable}."
    )

from xgboost import XGBClassifier
import matplotlib.pyplot as plt
from tqdm import tqdm

try:
    import cupy as cp
except ModuleNotFoundError:
    cp = None


def cuda_available() -> bool:
    if cp is None:
        return False
    try:
        return cp.cuda.runtime.getDeviceCount() > 0
    except Exception:
        return False


USE_CUDA = cuda_available()


def backend_asarray(a, dtype=None):
    if USE_CUDA:
        return cp.asarray(a, dtype=dtype)
    return np.asarray(a, dtype=dtype)


def backend_to_numpy(a):
    if USE_CUDA:
        return cp.asnumpy(a)
    return np.asarray(a)



# =========================
# 0) 配置
# =========================
ROOT_DATA5_CCDC2 = "/path/to/data2/CRC/CCDC2"
ML_DIR = os.path.join(ROOT_DATA5_CCDC2, "ML_results")

RANKED_FEATURES_TXT = os.path.join(ML_DIR, "global_ranked_features.txt")
MEAN_IMP_PATTERN = os.path.join(ML_DIR, "*_feature_mean_importance.csv")

# Evaluate the full K range by removing one feature at a time.
K_MIN = 1
K_STEP = 1   # 必须为 1 才是“逐个移除”

# XGBoost 参数（GPU 可用时走 GPU，否则自动回退 CPU）
XGB_PARAMS = dict(
    eval_metric="logloss",
    random_state=42,
    n_estimators=500,
    learning_rate=0.05,
    max_depth=4,
    subsample=0.8,
    colsample_bytree=0.8,
    reg_lambda=1.0,
    tree_method="hist",
    device="cuda:0" if USE_CUDA else "cpu",
    importance_type="gain",
    verbosity=0,
)

# Output paths.
OUT_FULL = os.path.join(ML_DIR, "rank_sum_curve_fullrange_xgb_gpu.csv")
OUT_AUC_BYK = os.path.join(ML_DIR, "auc_by_dataset_byK_fullrange_xgb_gpu.csv")
OUT_BEST_FEATURES = os.path.join(ML_DIR, "bestK_features_fullrange_xgb_gpu.txt")
OUT_SCATTER_PDF = os.path.join(ML_DIR, "paperstyle_GenomeNumber_vs_SumOfRank_xgb_gpu.pdf")


# =========================
# 0.1) 不同 dataset 的 seqkit_stats 根目录不一致：显式映射 + 自动探测
# =========================
SEQKIT_ROOT_BY_DATASET = {
    "FengQ_2015": "/path/to/data1",
    "YangJ_2020": "/path/to/data1/CRC_DATA",
    "YuJ_2015": "/path/to/data4/CRC_DATA",
    "YachidaS_2019": "/path/to/data4/CRC_DATA",
    "ONCOBIOME_2025": "/path/to/data1/CRC_DATA",
}
SEQKIT_ROOT_CANDIDATES = [
    "/path/to/data1/CRC_DATA",
    "/path/to/data4/CRC_DATA",
    "/path/to/data1",
]


# =========================
# 1) 工具函数
# =========================
def dataset_name_from_mean_imp(path: str) -> str:
    base = os.path.basename(path)
    suf = "_feature_mean_importance.csv"
    return base[:-len(suf)] if base.endswith(suf) else os.path.splitext(base)[0]


def load_ranked_features(path: str) -> list[str]:
    feats = []
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            s = line.strip()
            if s:
                feats.append(s)
    if not feats:
        raise ValueError(f"No features found in {path}")
    return feats


def resolve_seqkit_stats_file(dataset_name: str) -> str:
    tried = []

    if dataset_name in SEQKIT_ROOT_BY_DATASET:
        root = SEQKIT_ROOT_BY_DATASET[dataset_name]
        p = os.path.join(root, dataset_name, "seqkit_stats.txt")
        if os.path.exists(p):
            return p
        tried.append(p)

    for root in SEQKIT_ROOT_CANDIDATES:
        p = os.path.join(root, dataset_name, "seqkit_stats.txt")
        if os.path.exists(p):
            return p
        tried.append(p)

    raise FileNotFoundError(
        f"[ERROR] seqkit_stats.txt not found for dataset={dataset_name}\n"
        f"Tried:\n  - " + "\n  - ".join(tried)
    )


def resolve_group_file(dataset_name: str) -> str:
    candidates = [
        f"{dataset_name}_CRC_Group.txt",
        f"{dataset_name}_Group.txt",
        os.path.join(ML_DIR, f"{dataset_name}_CRC_Group.txt"),
        os.path.join(ML_DIR, f"{dataset_name}_Group.txt"),
    ]
    for p in candidates:
        if os.path.exists(p):
            return p
    raise FileNotFoundError(
        f"[ERROR] group file not found for dataset={dataset_name}\n"
        f"Tried:\n  - " + "\n  - ".join(candidates)
    )


def get_data_paths(dataset_name: str):
    data_dir = os.path.join(ROOT_DATA5_CCDC2, dataset_name, "DiTASiC", "abundance")
    seqkit_stats = resolve_seqkit_stats_file(dataset_name)
    group_file = resolve_group_file(dataset_name)
    return data_dir, seqkit_stats, group_file


def precompute_loo_splits_cp(n_samples: int, y_np: np.ndarray):
    """
    预先生成 LOOCV 的 (train_idx_cp, test_idx_cp, test_idx_int)
    这样每个 K / fold 不需要重复做 numpy<->cupy 转换。
    """
    loo = LeaveOneOut()
    dummy = np.empty((n_samples, 1), dtype=np.float32)
    splits = []
    for tr, te in loo.split(dummy, y_np):
        splits.append((backend_asarray(tr), backend_asarray(te), int(te[0])))
    return splits


# =========================
# 2) 加载 dataset 的 X/y，并对齐到 ranked_features（缺失补0）-> 转计算后端
# =========================
def load_dataset_Xy_aligned_gpu(dataset_name: str, ranked_features: list[str]):
    data_dir, seqkit_stats_file, group_file = get_data_paths(dataset_name)

    if not os.path.isdir(data_dir):
        raise FileNotFoundError(f"[ERROR] abundance dir not found: {data_dir}")

    # ---- group ----
    group_df = pd.read_csv(group_file, sep="\t", header=0)
    group_df.columns = ["Sample", "Group"]
    group_df["label"] = group_df["Group"].fillna("").astype(str).apply(lambda x: 1 if x.strip() == "CRC" else 0)
    group_df = group_df.set_index("Sample")

    # ---- seqkit ----
    seqkit_df = pd.read_csv(seqkit_stats_file, sep=r"\s+", header=0)
    seqkit_df["SampleID"] = seqkit_df["file"].str.split("/").str[0]
    seqkit_df["num_seqs"] = seqkit_df["num_seqs"].astype(str).str.replace(",", "").astype(int)
    seqkit_df["num_reads"] = seqkit_df["num_seqs"].astype(int)
    total_reads_per_sample = seqkit_df.groupby("SampleID")["num_reads"].sum()

    # ---- abundance ----
    file_list = glob(os.path.join(data_dir, "abundance_*.txt"))
    sample_data = []

    for file_path in file_list:
        sample_id = os.path.basename(file_path).replace("abundance_", "").replace(".txt", "")

        if sample_id not in group_df.index:
            continue
        if sample_id not in total_reads_per_sample.index:
            continue

        df = pd.read_csv(file_path, sep="\t")
        df = df[df["raw.pval"] <= 0.05]
        if df.empty:
            continue

        s = df.set_index("taxa.name")["count.estimate"]
        s.name = sample_id
        sample_data.append(s)

    if not sample_data:
        raise ValueError(f"{dataset_name}: no usable samples after filtering.")

    abundance_matrix = pd.DataFrame(sample_data).fillna(0)

    # reads 对齐
    reads_aligned = total_reads_per_sample.loc[abundance_matrix.index].astype(float)
    valid = reads_aligned > 0
    abundance_matrix = abundance_matrix.loc[valid]
    reads_aligned = reads_aligned.loc[valid]

    abundance_rel = abundance_matrix.div(reads_aligned, axis=0)

    # 对齐 label
    abundance_rel = abundance_rel.loc[group_df.index.intersection(abundance_rel.index)]
    y_np = group_df.loc[abundance_rel.index, "label"].to_numpy(dtype=int)

    # Align to global ranked_features, fill missing values with zero, and preserve column order.
    X_df = abundance_rel.reindex(columns=ranked_features, fill_value=0.0)
    X_np = X_df.to_numpy(dtype=np.float32)

    # ---- 转到当前计算后端（CUDA/CPU）----
    X_gpu = backend_asarray(X_np, dtype=np.float32)
    y_gpu = backend_asarray(y_np)

    # 预计算 LOOCV 索引
    splits_cp = precompute_loo_splits_cp(X_np.shape[0], y_np)

    return X_gpu, y_np, y_gpu, splits_cp


# =========================
# 3) 单个 dataset：LOOCV AUC（给定 K）
# =========================
def loocv_auc_xgb_gpu(X_gpu, y_np: np.ndarray, y_gpu, splits_cp, K: int) -> float:
    y_np = np.asarray(y_np).astype(int)
    if len(np.unique(y_np)) < 2:
        return np.nan

    oof = np.zeros(len(y_np), dtype=float)
    Xk = X_gpu[:, :K]  # 后端切片（保留前K个最重要特征）

    for tr_cp, te_cp, te_int in splits_cp:
        model = XGBClassifier(**XGB_PARAMS)
        model.fit(Xk[tr_cp], y_gpu[tr_cp])

        prob = backend_to_numpy(model.predict_proba(Xk[te_cp])[:, 1])[0]
        oof[te_int] = float(prob)

    return roc_auc_score(y_np, oof)


# =========================
# 4) Compute AUC and rank-sum across the full K range.
# =========================
def eval_rank_sum_fullrange_gpu(datasets_data: dict, K_list: list[int]):
    records = []

    # 外层：dataset 进度
    for ds_name, pack in tqdm(datasets_data.items(), total=len(datasets_data), desc="Datasets"):
        X_gpu, y_np, y_gpu, splits_cp = pack["X_gpu"], pack["y_np"], pack["y_gpu"], pack["splits_cp"]

        # 内层：K 进度（每个dataset一条）
        for K in tqdm(K_list, desc=f"K for {ds_name}", leave=False):
            auc_k = loocv_auc_xgb_gpu(X_gpu, y_np, y_gpu, splits_cp, K)
            records.append({"dataset": ds_name, "K": K, "AUC": auc_k})

    results = pd.DataFrame(records)

    ranked_parts = []
    for ds_name, sub in results.groupby("dataset"):
        sub = sub.copy()
        sub["AUC_rank_in_ds"] = rankdata(-sub["AUC"].values, method="average")
        ranked_parts.append(sub)
    ranked_df = pd.concat(ranked_parts, ignore_index=True)

    summary_df = (
        ranked_df.groupby("K")
        .agg(rank_sum=("AUC_rank_in_ds", "sum"), mean_auc=("AUC", "mean"))
        .reset_index()
    )

    best_row = summary_df.loc[summary_df["rank_sum"].idxmin()]
    best_K = int(best_row["K"])
    return best_K, summary_df, ranked_df



# =========================
# 5) 主流程
# =========================
def main():
    print(f"[INFO] backend: {'CUDA' if USE_CUDA else 'CPU'}")

    ranked_features = load_ranked_features(RANKED_FEATURES_TXT)
    print("Most important (first):", ranked_features[0])
    print("Least important (last):", ranked_features[-1])
    P = len(ranked_features)
    print(f"[INFO] ranked_features loaded: {P}")

    mean_imp_files = sorted(glob(MEAN_IMP_PATTERN))
    if not mean_imp_files:
        raise FileNotFoundError(f"No mean importance files found: {MEAN_IMP_PATTERN}")

    dataset_names = [dataset_name_from_mean_imp(p) for p in mean_imp_files]
    print("[INFO] datasets:", dataset_names)

    # 加载并转到当前计算后端
    datasets_data = {}
    for ds in dataset_names:
        X_gpu, y_np, y_gpu, splits_cp = load_dataset_Xy_aligned_gpu(ds, ranked_features)
        datasets_data[ds] = {"X_gpu": X_gpu, "y_np": y_np, "y_gpu": y_gpu, "splits_cp": splits_cp}
        print(f"[LOAD] {ds}: samples={X_gpu.shape[0]}, features_aligned={X_gpu.shape[1]} ({'CUDA' if USE_CUDA else 'CPU'})")

    # Remove features from least to most important by evaluating K from P to 1.
    K_list = list(range(P, K_MIN - 1, -K_STEP))
    print(f"[INFO] Evaluating full range K: {K_list[0]} -> {K_list[-1]} (step={K_STEP}), total={len(K_list)}")

    bestK, summary_df, ranked_df = eval_rank_sum_fullrange_gpu(datasets_data, K_list)

    # 保存结果
    ranked_df.to_csv(OUT_AUC_BYK, index=False)
    summary_df.to_csv(OUT_FULL, index=False)
    print(f"[SAVE] {OUT_AUC_BYK}")
    print(f"[SAVE] {OUT_FULL}")
    print(f"[INFO] bestK = {bestK}")

    # 输出 bestK features（保留前 bestK 个最重要特征）
    best_features = ranked_features[:bestK]
    with open(OUT_BEST_FEATURES, "w", encoding="utf-8") as f:
        for feat in best_features:
            f.write(feat + "\n")
    print(f"[SAVE] bestK features list -> {OUT_BEST_FEATURES}")

    # =========================
    #  - 横轴从 0 开始到 P
    #  - 黑色小点
    #  - y 轴是 Sum of rank
    # =========================
    plot_df = summary_df.sort_values("K")  # Sort K in ascending order for the left-to-right axis.

    plt.rcParams.update({
        "font.size": 14,
        "axes.linewidth": 1.0,
    })

    plt.figure(figsize=(6.0, 4.5))
    plt.scatter(plot_df["K"], plot_df["rank_sum"], s=8, c="black", marker=".", linewidths=0)

    plt.xlabel("Genome Number")   
    plt.ylabel("Sum of rank")    
    plt.xlim(0, P)               
    plt.grid(False)

    plt.tight_layout()
    plt.savefig(OUT_SCATTER_PDF)
    plt.close()
    print(f"[SAVE] paper-style scatter -> {OUT_SCATTER_PDF}")

    print("\n[TOP 10 by smallest rank_sum]:")
    print(summary_df.sort_values("rank_sum").head(10))


if __name__ == "__main__":
    main()
