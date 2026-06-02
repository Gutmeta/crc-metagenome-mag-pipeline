import os
from glob import glob
import pandas as pd
import numpy as np

ML_DIR = "/path/to/data2/CRC/CCDC2/ML_results"
PATTERN = "*_feature_mean_importance.csv"

OUT_GLOBAL = os.path.join(ML_DIR, "global_feature_rank_from_mean_importance.csv")
OUT_RANKED_LIST = os.path.join(ML_DIR, "global_ranked_features.txt")

def dataset_name_from_file(path: str) -> str:
    base = os.path.basename(path)
    suf = "_feature_mean_importance.csv"
    return base[:-len(suf)] if base.endswith(suf) else os.path.splitext(base)[0]

def load_mean_importance(path: str) -> pd.Series:
    df = pd.read_csv(path)
    df = df[["Feature", "MeanImportance"]].copy()
    df["MeanImportance"] = pd.to_numeric(df["MeanImportance"], errors="coerce")
    df = df.dropna(subset=["Feature", "MeanImportance"])
    # 若有重复特征名，取均值
    return df.groupby("Feature")["MeanImportance"].mean()

def infer_weight_from_prediction_csv(ds_name: str) -> int:
    pred_path = os.path.join(ML_DIR, f"{ds_name}_prediction_results.csv")
    if os.path.exists(pred_path):
        try:
            return int(pd.read_csv(pred_path).shape[0])
        except Exception:
            pass
    # 找不到就退化为 1（也能跑）
    return 1

files = sorted(glob(os.path.join(ML_DIR, PATTERN)))
if not files:
    raise FileNotFoundError(f"No files matched {os.path.join(ML_DIR, PATTERN)}")

# 读取每个dataset的 mean importance + 权重
imp = {}
weights = {}
for f in files:
    ds = dataset_name_from_file(f)
    imp[ds] = load_mean_importance(f)
    weights[ds] = infer_weight_from_prediction_csv(ds)

# union 对齐（缺失补0）
all_features = sorted(set().union(*[set(s.index) for s in imp.values()]))
imp_df = pd.DataFrame({ds: s.reindex(all_features, fill_value=0.0) for ds, s in imp.items()})
w = pd.Series(weights, dtype=float)

# 1) 论文式：across all models 的加权平均重要性（权重≈fold数≈样本数）
global_mean_imp = (imp_df * w).sum(axis=1) / w.sum()

# 2) 平均名次（更鲁棒）：每个dataset内部按 importance 排名，然后按权重平均
rank_df = imp_df.rank(axis=0, ascending=False, method="average")  # 1=最重要
global_avg_rank = (rank_df * w).sum(axis=1) / w.sum()

out = pd.DataFrame({
    "Feature": all_features,
    "GlobalMeanImportance_weighted": global_mean_imp.values,
    "GlobalAvgRank_weighted": global_avg_rank.values,
}).set_index("Feature")

# 给两套全局排名
out["GlobalRank_ByMeanImportance"] = out["GlobalMeanImportance_weighted"].rank(ascending=False, method="min").astype(int)
out["GlobalRank_ByAvgRank"] = out["GlobalAvgRank_weighted"].rank(ascending=True, method="min").astype(int)

# 可选：把每个dataset的 mean importance 也拼进去便于检查
out = out.join(imp_df.add_prefix("MeanImp__"), how="left")

# 排序：优先用平均名次（更鲁棒），再用平均重要性
out = out.sort_values(["GlobalRank_ByAvgRank", "GlobalRank_ByMeanImportance"])

out.to_csv(OUT_GLOBAL)

ranked_features = out.index.tolist()
with open(OUT_RANKED_LIST, "w", encoding="utf-8") as f:
    for feat in ranked_features:
        f.write(feat + "\n")

print("[OK] Saved:", OUT_GLOBAL)
print("[OK] Saved:", OUT_RANKED_LIST)
print("Datasets loaded:", list(imp.keys()))
print("Weights (approx folds):", weights)
print("Top 20 features:", ranked_features[:20])
