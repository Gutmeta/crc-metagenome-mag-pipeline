#!/usr/bin/env python3
"""Unified leave-one-dataset-out best-k AUC analysis for CRA, CRC, and IBD."""

from __future__ import annotations

import math
import os
from dataclasses import dataclass
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", str(Path(__file__).resolve().parent / ".mplconfig"))
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.preprocessing import StandardScaler


ROOT = Path(__file__).resolve().parents[1]
OUTDIR = ROOT / "requested_topjournal_figures"
INVENTORY = ROOT / "key_species_overlap_v2" / "disease_arm_inventory.tsv"

DISEASES = ("CRA", "CRC", "IBD")
CONTROL_STATUSES = {"ctr", "control patient"}
IBD_CASE_STATUSES = {
    "crohn's disease",
    "ulcerative colitis",
    "inflammatory bowel disease",
    "indeterminate colitis",
}

MIN_CANDIDATE_K = 1
MAX_CANDIDATE_K = 200
MIN_DETECTION_FRAC = 0.05
MIN_DETECTION_N = 10
MIN_MEAN_ABUND = 1e-5
LOG_PSEUDOCOUNT = 1e-6
CORR_CUTOFF = 0.90
MAX_RANK_POOL = 1200
N_AUC_BOOTSTRAP = 5000
N_AUC_PERMUTATION = 5000
RANDOM_SEED = 20260524


@dataclass
class ArmData:
    disease: str
    study_code: str
    arm_id: str
    metadata_path: Path
    profile_path: Path
    x_raw: pd.DataFrame
    y: pd.Series

    @property
    def n_case(self) -> int:
        return int(self.y.sum())

    @property
    def n_control(self) -> int:
        return int((self.y == 0).sum())


def label_status(disease: str, status: str) -> str | None:
    s = str(status or "").strip().lower()
    if s in CONTROL_STATUSES:
        return "control"
    if disease == "CRA":
        if "adenoma" in s or "adenomatous" in s:
            return "case"
        return None
    if disease == "CRC":
        return "case" if s == "colorectal cancer" else None
    if disease == "IBD":
        return "case" if s in IBD_CASE_STATUSES else None
    return None


def safe_float(value: object, default: float = np.nan) -> float:
    try:
        return float(value)
    except Exception:
        return default


def is_true(value: object) -> bool:
    return str(value).strip().lower() == "true"


def resolve_input_path(value: object) -> Path:
    path = Path(str(value))
    if path.exists():
        return path
    if path.is_absolute():
        for anchor in ("CRA", "CRC", "IBD", "IBS"):
            if anchor in path.parts:
                candidate = ROOT.joinpath(*path.parts[path.parts.index(anchor) :])
                if candidate.exists():
                    return candidate
    return path


def load_metadata(path: Path, disease: str) -> pd.DataFrame:
    meta = pd.read_csv(path, sep="\t", dtype=str)
    required = {"sample_alias", "subject_disease_status"}
    missing = required.difference(meta.columns)
    if missing:
        raise ValueError(f"{path} missing columns: {sorted(missing)}")
    meta = meta[["sample_alias", "subject_disease_status"]].copy()
    meta["group"] = meta["subject_disease_status"].map(lambda x: label_status(disease, x))
    meta = meta.dropna(subset=["group"])
    meta = meta.drop_duplicates("sample_alias", keep="first")
    meta["y"] = (meta["group"] == "case").astype(int)
    return meta[["sample_alias", "y"]]


def load_metaphlan_species(path: Path) -> pd.DataFrame:
    df = pd.read_csv(
        path,
        sep="\t",
        compression="gzip",
        usecols=["sample_alias", "clade_name", "rel_abund"],
        dtype={"sample_alias": str, "clade_name": str, "rel_abund": np.float32},
    )
    df = df[df["clade_name"].str.startswith("s__", na=False)]
    mat = df.pivot_table(
        index="sample_alias",
        columns="clade_name",
        values="rel_abund",
        aggfunc="sum",
        fill_value=0.0,
    )
    return mat.astype(np.float32)


def load_arm(row: pd.Series) -> ArmData:
    disease = str(row["disease"])
    study_code = str(row["study_code"])
    arm_id = str(row["arm_id"])
    metadata_path = resolve_input_path(row["metadata_path"])
    profile_path = resolve_input_path(row["profile_path"])

    meta = load_metadata(metadata_path, disease)
    x_raw = load_metaphlan_species(profile_path)
    common = meta["sample_alias"][meta["sample_alias"].isin(x_raw.index)]
    meta = meta.set_index("sample_alias").loc[common]
    x_raw = x_raw.loc[common]
    y = meta["y"].astype(int)
    if y.nunique() != 2:
        raise ValueError(f"{arm_id} has only one class after profile matching")
    return ArmData(disease, study_code, arm_id, metadata_path, profile_path, x_raw, y)


def combine_arms(arms: list[ArmData]) -> tuple[pd.DataFrame, pd.Series, pd.Series]:
    features = sorted(set().union(*(set(arm.x_raw.columns) for arm in arms)))
    frames = []
    ys = []
    cohorts = []
    for arm in arms:
        frame = arm.x_raw.reindex(columns=features, fill_value=0.0)
        frame.index = [f"{arm.study_code}|{sample}" for sample in frame.index]
        frames.append(frame)
        y = arm.y.copy()
        y.index = frame.index
        ys.append(y)
        cohorts.extend([arm.study_code] * frame.shape[0])
    x = pd.concat(frames, axis=0)
    y = pd.concat(ys, axis=0).astype(int)
    cohort = pd.Series(cohorts, index=x.index, name="cohort")
    return x, y, cohort


def mannwhitney_p(case_values: np.ndarray, control_values: np.ndarray) -> float:
    try:
        if np.all(case_values == case_values[0]) and np.all(control_values == control_values[0]):
            if case_values[0] == control_values[0]:
                return 1.0
        return float(mannwhitneyu(case_values, control_values, alternative="two-sided").pvalue)
    except Exception:
        return 1.0


def rank_species(x_train_raw: pd.DataFrame, y_train: pd.Series) -> pd.DataFrame:
    detected = (x_train_raw > 0).sum(axis=0)
    prevalence = detected / max(1, x_train_raw.shape[0])
    mean_abund = x_train_raw.mean(axis=0)
    keep = ((prevalence >= MIN_DETECTION_FRAC) | (detected >= MIN_DETECTION_N)) & (mean_abund >= MIN_MEAN_ABUND)
    kept = x_train_raw.loc[:, keep]
    if kept.shape[1] == 0:
        return pd.DataFrame(columns=["species", "score", "log2fc", "p_value", "detected", "mean_abund"])

    y_values = y_train.to_numpy()
    case_mask = y_values == 1
    control_mask = y_values == 0
    case_mean = kept.iloc[case_mask].mean(axis=0).to_numpy()
    control_mean = kept.iloc[control_mask].mean(axis=0).to_numpy()
    log2fc = np.log2((case_mean + LOG_PSEUDOCOUNT) / (control_mean + LOG_PSEUDOCOUNT))

    x_log = np.log10(kept.to_numpy(dtype=np.float32) + LOG_PSEUDOCOUNT)
    p_values = np.array(
        [mannwhitney_p(x_log[case_mask, i], x_log[control_mask, i]) for i in range(x_log.shape[1])],
        dtype=float,
    )
    p_values = np.where(np.isfinite(p_values), np.clip(p_values, 1e-300, 1.0), 1.0)
    score = -np.log10(p_values) * (np.abs(log2fc) + 1e-6)

    ranked = pd.DataFrame(
        {
            "species": kept.columns,
            "score": score,
            "log2fc": log2fc,
            "p_value": p_values,
            "detected": detected.loc[kept.columns].to_numpy(),
            "mean_abund": mean_abund.loc[kept.columns].to_numpy(),
        }
    )
    return ranked.sort_values(["score", "p_value", "species"], ascending=[False, True, True]).reset_index(drop=True)


def corr_deduplicate(
    x_train_raw: pd.DataFrame,
    ranked: pd.DataFrame,
    max_k: int,
    corr_cutoff: float = CORR_CUTOFF,
) -> pd.DataFrame:
    if ranked.empty:
        return ranked.copy()
    pool = ranked.head(max(MAX_RANK_POOL, max_k)).copy()
    candidates = pool["species"].tolist()
    x_log = np.log10(x_train_raw.reindex(columns=candidates, fill_value=0.0).to_numpy(dtype=np.float32) + LOG_PSEUDOCOUNT)
    selected_indices: list[int] = []
    selected_species: list[str] = []
    for idx, species in enumerate(candidates):
        if len(selected_species) >= max_k:
            break
        if not selected_indices:
            selected_indices.append(idx)
            selected_species.append(species)
            continue
        values = x_log[:, idx]
        selected_values = x_log[:, selected_indices]
        value_sd = values.std()
        selected_sd = selected_values.std(axis=0)
        valid = (value_sd > 0) & (selected_sd > 0)
        if valid.any():
            corrs = np.corrcoef(values, selected_values[:, valid], rowvar=False)[0, 1:]
            max_corr = np.nanmax(np.abs(corrs)) if corrs.size else 0.0
        else:
            max_corr = 0.0
        if not np.isfinite(max_corr) or max_corr < corr_cutoff:
            selected_indices.append(idx)
            selected_species.append(species)
    selected = ranked[ranked["species"].isin(selected_species)].copy()
    selected["feature_rank"] = selected["species"].map({sp: i + 1 for i, sp in enumerate(selected_species)})
    return selected.sort_values("feature_rank").reset_index(drop=True)


def fit_predict_scores(
    x_train_raw: pd.DataFrame,
    y_train: pd.Series,
    x_test_raw: pd.DataFrame,
    y_test: pd.Series,
    features: list[str],
) -> tuple[float, np.ndarray]:
    if not features or y_test.nunique() != 2 or y_train.nunique() != 2:
        return np.nan, np.full(y_test.shape[0], np.nan, dtype=float)
    x_train = np.log10(x_train_raw.reindex(columns=features, fill_value=0.0).to_numpy(dtype=np.float32) + LOG_PSEUDOCOUNT)
    x_test = np.log10(x_test_raw.reindex(columns=features, fill_value=0.0).to_numpy(dtype=np.float32) + LOG_PSEUDOCOUNT)
    scaler = StandardScaler()
    x_train = scaler.fit_transform(x_train)
    x_test = scaler.transform(x_test)
    model = LogisticRegression(
        class_weight="balanced",
        solver="liblinear",
        C=1.0,
        max_iter=1000,
        random_state=RANDOM_SEED,
    )
    model.fit(x_train, y_train.to_numpy())
    pred = model.predict_proba(x_test)[:, 1]
    return float(roc_auc_score(y_test.to_numpy(), pred)), pred


def fit_predict_auc(
    x_train_raw: pd.DataFrame,
    y_train: pd.Series,
    x_test_raw: pd.DataFrame,
    y_test: pd.Series,
    features: list[str],
) -> float:
    auc, _ = fit_predict_scores(x_train_raw, y_train, x_test_raw, y_test, features)
    return auc


def bootstrap_auc_ci(
    y_true: pd.Series,
    pred: np.ndarray,
    n_bootstrap: int = N_AUC_BOOTSTRAP,
    seed: int = RANDOM_SEED,
) -> tuple[float, float]:
    y = y_true.to_numpy(dtype=int)
    pred = np.asarray(pred, dtype=float)
    case_idx = np.where(y == 1)[0]
    control_idx = np.where(y == 0)[0]
    if len(case_idx) == 0 or len(control_idx) == 0:
        return np.nan, np.nan
    rng = np.random.default_rng(seed)
    aucs = np.empty(n_bootstrap, dtype=float)
    for i in range(n_bootstrap):
        sampled_case = rng.choice(case_idx, size=len(case_idx), replace=True)
        sampled_control = rng.choice(control_idx, size=len(control_idx), replace=True)
        idx = np.concatenate([sampled_case, sampled_control])
        aucs[i] = roc_auc_score(y[idx], pred[idx])
    return float(np.nanpercentile(aucs, 2.5)), float(np.nanpercentile(aucs, 97.5))


def permutation_auc_p(
    y_true: pd.Series,
    pred: np.ndarray,
    observed_auc: float,
    n_permutation: int = N_AUC_PERMUTATION,
    seed: int = RANDOM_SEED,
) -> float:
    y = y_true.to_numpy(dtype=int)
    pred = np.asarray(pred, dtype=float)
    if not np.isfinite(observed_auc) or len(np.unique(y)) != 2:
        return np.nan
    rng = np.random.default_rng(seed)
    hits = 0
    for _ in range(n_permutation):
        shuffled = rng.permutation(y)
        perm_auc = roc_auc_score(shuffled, pred)
        if perm_auc >= observed_auc:
            hits += 1
    return float((hits + 1) / (n_permutation + 1))


def disease_bootstrap_summary(source: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    rng = np.random.default_rng(RANDOM_SEED)
    for disease, group in source.groupby("disease", sort=False):
        aucs = group["auc"].to_numpy(dtype=float)
        if aucs.size == 0:
            continue
        boot_median = np.empty(N_AUC_BOOTSTRAP, dtype=float)
        boot_mean = np.empty(N_AUC_BOOTSTRAP, dtype=float)
        for i in range(N_AUC_BOOTSTRAP):
            sample = rng.choice(aucs, size=aucs.size, replace=True)
            boot_median[i] = np.median(sample)
            boot_mean[i] = np.mean(sample)
        rows.append(
            {
                "disease": disease,
                "n_cohorts": int(aucs.size),
                "median_auc": float(np.median(aucs)),
                "median_auc_ci_low": float(np.percentile(boot_median, 2.5)),
                "median_auc_ci_high": float(np.percentile(boot_median, 97.5)),
                "mean_auc": float(np.mean(aucs)),
                "mean_auc_ci_low": float(np.percentile(boot_mean, 2.5)),
                "mean_auc_ci_high": float(np.percentile(boot_mean, 97.5)),
                "min_auc": float(np.min(aucs)),
                "max_auc": float(np.max(aucs)),
                "median_best_k": float(np.median(group["best_k"].to_numpy(dtype=float))),
            }
        )
    return pd.DataFrame(rows)


def choose_final_k_from_curves(curve_all: pd.DataFrame, disease: str, fallback_k: int) -> tuple[int, pd.DataFrame]:
    disease_curve = curve_all[curve_all["disease"] == disease].copy()
    if disease_curve.empty:
        summary = pd.DataFrame(
            [
                {
                    "disease": disease,
                    "candidate_k": int(fallback_k),
                    "mean_inner_auc": np.nan,
                    "sd_inner_auc": np.nan,
                    "n_inner_evaluations": 0,
                    "selected_final_k": True,
                }
            ]
        )
        return int(fallback_k), summary
    summary = (
        disease_curve.groupby("candidate_k", as_index=False)
        .agg(
            mean_inner_auc=("inner_auc", "mean"),
            sd_inner_auc=("inner_auc", "std"),
            n_inner_evaluations=("inner_auc", "count"),
        )
        .sort_values(["mean_inner_auc", "candidate_k"], ascending=[False, True])
        .reset_index(drop=True)
    )
    final_k = int(summary.iloc[0]["candidate_k"])
    summary.insert(0, "disease", disease)
    summary["selected_final_k"] = summary["candidate_k"] == final_k
    return final_k, summary


def candidate_k_values(max_effective_k: int = MAX_CANDIDATE_K) -> range:
    max_k = max(MIN_CANDIDATE_K, min(MAX_CANDIDATE_K, int(max_effective_k)))
    return range(MIN_CANDIDATE_K, max_k + 1)


def nested_selection_stats(feature_rows: list[dict[str, object]], source: pd.DataFrame) -> pd.DataFrame:
    feature_df = pd.DataFrame(feature_rows)
    if feature_df.empty:
        return pd.DataFrame(columns=["disease", "species", "selection_count", "selection_frequency", "median_nested_feature_rank"])
    n_folds = source[source["status"] == "ok"].groupby("disease")["left_out_cohort"].nunique().to_dict()
    stats = (
        feature_df.groupby(["disease", "species"], as_index=False)
        .agg(selection_count=("left_out_cohort", "nunique"), median_nested_feature_rank=("feature_rank", "median"))
    )
    stats["selection_frequency"] = stats.apply(
        lambda row: float(row["selection_count"]) / max(1, int(n_folds.get(row["disease"], 0))),
        axis=1,
    )
    return stats


def compute_final_signatures(
    arms_by_disease: dict[str, list[ArmData]],
    source: pd.DataFrame,
    feature_rows: list[dict[str, object]],
    curve_all: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    selection_stats = nested_selection_stats(feature_rows, source)
    signature_feature_rows: list[dict[str, object]] = []
    signature_auc_rows: list[dict[str, object]] = []
    signature_prediction_rows: list[dict[str, object]] = []
    k_summary_rows: list[pd.DataFrame] = []

    for disease in DISEASES:
        arms = arms_by_disease[disease]
        if len(arms) < 3:
            continue
        fallback_k = int(source.loc[source["disease"] == disease, "best_k"].median())
        final_k, k_summary = choose_final_k_from_curves(curve_all, disease, fallback_k)
        k_summary_rows.append(k_summary)

        x_full, y_full, _ = combine_arms(arms)
        ranked_full = rank_species(x_full, y_full)
        selected_full = corr_deduplicate(x_full, ranked_full, max_k=MAX_CANDIDATE_K)
        final_panel = selected_full.head(final_k).copy()
        features = final_panel["species"].tolist()
        final_effective_k = len(features)
        if final_effective_k == 0:
            continue

        stats_d = selection_stats[selection_stats["disease"] == disease].set_index("species")
        k_selected = k_summary[k_summary["selected_final_k"]].iloc[0]
        for _, feat in final_panel.iterrows():
            sp = feat["species"]
            signature_feature_rows.append(
                {
                    "disease": disease,
                    "final_k": final_k,
                    "final_effective_k": final_effective_k,
                    "feature_rank": int(feat["feature_rank"]),
                    "species": sp,
                    "full_disease_score": safe_float(feat["score"]),
                    "full_disease_log2fc": safe_float(feat["log2fc"]),
                    "full_disease_p_value": safe_float(feat["p_value"]),
                    "full_disease_detected": int(feat["detected"]),
                    "full_disease_mean_abund": safe_float(feat["mean_abund"]),
                    "nested_selection_count": int(stats_d.loc[sp, "selection_count"]) if sp in stats_d.index else 0,
                    "nested_selection_frequency": safe_float(stats_d.loc[sp, "selection_frequency"], 0.0) if sp in stats_d.index else 0.0,
                    "median_nested_feature_rank": safe_float(stats_d.loc[sp, "median_nested_feature_rank"]) if sp in stats_d.index else np.nan,
                    "final_k_mean_inner_auc": safe_float(k_selected["mean_inner_auc"]),
                    "final_k_sd_inner_auc": safe_float(k_selected["sd_inner_auc"]),
                    "final_k_n_inner_evaluations": int(k_selected["n_inner_evaluations"]),
                    "feature_selection_scope": "full_disease_post_selection",
                    "model_family": "balanced_logistic_regression",
                }
            )

        print(
            f"[{disease}] locked final signature: final_k={final_k}, effective_k={final_effective_k}, "
            f"mean inner AUC={safe_float(k_selected['mean_inner_auc']):.3f}",
            flush=True,
        )

        for left_out in arms:
            train_arms = [arm for arm in arms if arm.arm_id != left_out.arm_id]
            x_train, y_train, _ = combine_arms(train_arms)
            x_test, y_test, _ = combine_arms([left_out])
            train_case, train_control = summarize_classes(y_train)
            test_case, test_control = summarize_classes(y_test)
            auc, pred = fit_predict_scores(x_train, y_train, x_test, y_test, features)
            row_seed = RANDOM_SEED + 1009 * len(signature_auc_rows) + 503
            ci_low, ci_high = bootstrap_auc_ci(y_test, pred, seed=row_seed)
            perm_p = permutation_auc_p(y_test, pred, auc, seed=row_seed + 37)
            signature_auc_rows.append(
                {
                    "disease": disease,
                    "left_out_cohort": left_out.study_code,
                    "arm_id": left_out.arm_id,
                    "status": "ok" if np.isfinite(auc) else "failed_auc",
                    "auc": auc,
                    "auc_ci_low": ci_low,
                    "auc_ci_high": ci_high,
                    "auc_bootstrap_n": N_AUC_BOOTSTRAP,
                    "auc_permutation_p_one_sided": perm_p,
                    "auc_permutation_n": N_AUC_PERMUTATION,
                    "significant_auc_gt_0_5": bool(np.isfinite(perm_p) and perm_p < 0.05 and auc > 0.5),
                    "best_k": final_k,
                    "final_k": final_k,
                    "effective_k": final_effective_k,
                    "n_train": int(y_train.shape[0]),
                    "n_test": int(y_test.shape[0]),
                    "train_case": train_case,
                    "train_control": train_control,
                    "test_case": test_case,
                    "test_control": test_control,
                    "selected_features": ";".join(features),
                    "feature_selection_scope": "full_disease_post_selection",
                    "model_training_scope": "leave_one_cohort_out_after_locked_panel",
                    "held_out_cohort_excluded_from_model_training": True,
                    "classifier": "balanced_logistic_regression",
                    "transform": "log10(relative_abundance + 1e-6), training-fold standard scaling",
                }
            )
            for sample_index, y_value, pred_value in zip(x_test.index, y_test.to_numpy(dtype=int), pred):
                sample_name = str(sample_index).split("|", 1)[1] if "|" in str(sample_index) else str(sample_index)
                signature_prediction_rows.append(
                    {
                        "disease": disease,
                        "left_out_cohort": left_out.study_code,
                        "sample_alias": sample_name,
                        "y_true": int(y_value),
                        "predicted_probability": float(pred_value),
                        "final_k": final_k,
                        "effective_k": final_effective_k,
                        "feature_selection_scope": "full_disease_post_selection",
                    }
                )
            print(
                f"  locked panel left out {left_out.study_code}: AUC={auc:.3f}, final_k={final_k}",
                flush=True,
            )

    feature_table = pd.DataFrame(signature_feature_rows)
    auc_table = pd.DataFrame(signature_auc_rows)
    prediction_table = pd.DataFrame(signature_prediction_rows)
    k_summary_table = pd.concat(k_summary_rows, ignore_index=True) if k_summary_rows else pd.DataFrame()
    disease_summary = disease_bootstrap_summary(auc_table) if not auc_table.empty else pd.DataFrame()
    if not disease_summary.empty:
        final_k_map = auc_table.groupby("disease")["final_k"].first().to_dict()
        effective_k_map = auc_table.groupby("disease")["effective_k"].first().to_dict()
        disease_summary["final_k"] = disease_summary["disease"].map(final_k_map).astype(int)
        disease_summary["final_effective_k"] = disease_summary["disease"].map(effective_k_map).astype(int)
        disease_summary["feature_selection_scope"] = "full_disease_post_selection"
    return feature_table, auc_table, prediction_table, disease_summary, k_summary_table


def choose_best_k(train_arms: list[ArmData]) -> tuple[int, pd.DataFrame]:
    rows: list[dict[str, object]] = []
    for validation_arm in train_arms:
        inner_train = [arm for arm in train_arms if arm.arm_id != validation_arm.arm_id]
        if len(inner_train) < 2:
            continue
        x_inner_train, y_inner_train, _ = combine_arms(inner_train)
        x_val, y_val, _ = combine_arms([validation_arm])
        ranked = rank_species(x_inner_train, y_inner_train)
        selected = corr_deduplicate(x_inner_train, ranked, max_k=MAX_CANDIDATE_K)
        selected_features = selected["species"].tolist()
        for k in candidate_k_values(len(selected_features)):
            features = selected_features[:k]
            auc = fit_predict_auc(x_inner_train, y_inner_train, x_val, y_val, features)
            rows.append(
                {
                    "inner_validation_cohort": validation_arm.study_code,
                    "candidate_k": k,
                    "effective_k": len(features),
                    "inner_auc": auc,
                    "n_inner_train": int(y_inner_train.shape[0]),
                    "n_inner_validation": int(y_val.shape[0]),
                }
            )
    curve = pd.DataFrame(rows)
    if curve.empty:
        return MIN_CANDIDATE_K, curve
    summary = (
        curve.groupby("candidate_k", as_index=False)
        .agg(mean_inner_auc=("inner_auc", "mean"), sd_inner_auc=("inner_auc", "std"), n_inner_folds=("inner_auc", "count"))
        .sort_values(["mean_inner_auc", "candidate_k"], ascending=[False, True])
    )
    best_k = int(summary.iloc[0]["candidate_k"])
    curve = curve.merge(summary, on="candidate_k", how="left")
    curve["selected_best_k"] = curve["candidate_k"] == best_k
    return best_k, curve


def summarize_classes(y: pd.Series) -> tuple[int, int]:
    return int(y.sum()), int((y == 0).sum())


def run() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    inventory = pd.read_csv(INVENTORY, sep="\t")
    inventory = inventory[
        inventory["disease"].isin(DISEASES)
        & inventory["included_main"].map(is_true)
        & inventory["profile_path"].notna()
        & (inventory["profile_path"].astype(str) != "")
    ].copy()

    arms_by_disease: dict[str, list[ArmData]] = {disease: [] for disease in DISEASES}
    exclusions: list[dict[str, object]] = []
    for _, row in inventory.iterrows():
        try:
            arm = load_arm(row)
            arms_by_disease[arm.disease].append(arm)
        except Exception as exc:
            exclusions.append(
                {
                    "disease": row.get("disease", ""),
                    "study_code": row.get("study_code", ""),
                    "arm_id": row.get("arm_id", ""),
                    "status": "failed_loading",
                    "message": str(exc),
                }
            )

    source_rows: list[dict[str, object]] = []
    feature_rows: list[dict[str, object]] = []
    prediction_rows: list[dict[str, object]] = []
    curve_rows: list[pd.DataFrame] = []

    for disease in DISEASES:
        arms = arms_by_disease[disease]
        if len(arms) < 3:
            for arm in arms:
                exclusions.append(
                    {
                        "disease": disease,
                        "study_code": arm.study_code,
                        "arm_id": arm.arm_id,
                        "status": "excluded",
                        "message": f"fewer_than_3_eligible_arms_for_lodo(n={len(arms)})",
                    }
                )
            continue
        print(f"[{disease}] running {len(arms)} outer LODO folds", flush=True)
        for left_out in arms:
            train_arms = [arm for arm in arms if arm.arm_id != left_out.arm_id]
            x_train, y_train, _ = combine_arms(train_arms)
            x_test, y_test, _ = combine_arms([left_out])
            train_case, train_control = summarize_classes(y_train)
            test_case, test_control = summarize_classes(y_test)
            if min(train_case, train_control, test_case, test_control) < 1:
                exclusions.append(
                    {
                        "disease": disease,
                        "study_code": left_out.study_code,
                        "arm_id": left_out.arm_id,
                        "status": "excluded",
                        "message": "one_or_more_lodo_sets_have_single_class",
                    }
                )
                continue

            best_k, curve = choose_best_k(train_arms)
            if not curve.empty:
                curve.insert(0, "left_out_cohort", left_out.study_code)
                curve.insert(0, "disease", disease)
                curve_rows.append(curve)

            ranked = rank_species(x_train, y_train)
            selected = corr_deduplicate(x_train, ranked, max_k=MAX_CANDIDATE_K)
            features = selected["species"].head(best_k).tolist()
            auc, pred = fit_predict_scores(x_train, y_train, x_test, y_test, features)
            ci_low, ci_high = bootstrap_auc_ci(
                y_test,
                pred,
                seed=RANDOM_SEED + 101 * len(source_rows) + 7,
            )
            perm_p = permutation_auc_p(
                y_test,
                pred,
                auc,
                seed=RANDOM_SEED + 101 * len(source_rows) + 31,
            )
            status = "ok" if np.isfinite(auc) else "failed_auc"
            source_rows.append(
                {
                    "disease": disease,
                    "left_out_cohort": left_out.study_code,
                    "arm_id": left_out.arm_id,
                    "status": status,
                    "auc": auc,
                    "auc_ci_low": ci_low,
                    "auc_ci_high": ci_high,
                    "auc_bootstrap_n": N_AUC_BOOTSTRAP,
                    "auc_permutation_p_one_sided": perm_p,
                    "auc_permutation_n": N_AUC_PERMUTATION,
                    "significant_auc_gt_0_5": bool(np.isfinite(perm_p) and perm_p < 0.05 and auc > 0.5),
                    "best_k": best_k,
                    "effective_k": len(features),
                    "n_train": int(y_train.shape[0]),
                    "n_test": int(y_test.shape[0]),
                    "train_case": train_case,
                    "train_control": train_control,
                    "test_case": test_case,
                    "test_control": test_control,
                    "n_train_species_before_filter": int(x_train.shape[1]),
                    "n_ranked_after_filter": int(ranked.shape[0]),
                    "n_after_correlation_filter": int(selected.shape[0]),
                    "selected_features": ";".join(features),
                    "classifier": "balanced_logistic_regression",
                    "transform": "log10(relative_abundance + 1e-6), training-fold standard scaling",
                }
            )
            for sample_index, y_value, pred_value in zip(x_test.index, y_test.to_numpy(dtype=int), pred):
                sample_name = str(sample_index).split("|", 1)[1] if "|" in str(sample_index) else str(sample_index)
                prediction_rows.append(
                    {
                        "disease": disease,
                        "left_out_cohort": left_out.study_code,
                        "sample_alias": sample_name,
                        "y_true": int(y_value),
                        "predicted_probability": float(pred_value),
                        "best_k": best_k,
                    }
                )
            for _, feat in selected.head(best_k).iterrows():
                feature_rows.append(
                    {
                        "disease": disease,
                        "left_out_cohort": left_out.study_code,
                        "feature_rank": int(feat["feature_rank"]),
                        "species": feat["species"],
                        "score": safe_float(feat["score"]),
                        "log2fc": safe_float(feat["log2fc"]),
                        "p_value": safe_float(feat["p_value"]),
                        "detected": int(feat["detected"]),
                        "mean_abund": safe_float(feat["mean_abund"]),
                    }
                )
            print(
                f"  left out {left_out.study_code}: AUC={auc:.3f}, best_k={best_k}, effective_k={len(features)}",
                flush=True,
            )

    source = pd.DataFrame(source_rows)
    if not source.empty:
        source["auc"] = source["auc"].astype(float)
    feature_df = pd.DataFrame(feature_rows)
    prediction_df = pd.DataFrame(prediction_rows)
    source.to_csv(OUTDIR / "figure1_source_data.tsv", sep="\t", index=False)
    feature_df.to_csv(OUTDIR / "figure1_selected_features.tsv", sep="\t", index=False)
    prediction_df.to_csv(OUTDIR / "figure1_predictions.tsv", sep="\t", index=False)
    disease_bootstrap_summary(source).to_csv(OUTDIR / "figure1_disease_summary.tsv", sep="\t", index=False)
    if curve_rows:
        curve_all = pd.concat(curve_rows, ignore_index=True)
        curve_all.to_csv(OUTDIR / "figure1_bestk_curve.tsv", sep="\t", index=False)
    else:
        curve_all = pd.DataFrame()
        pd.DataFrame().to_csv(OUTDIR / "figure1_bestk_curve.tsv", sep="\t", index=False)

    (
        final_features,
        final_auc,
        final_predictions,
        final_disease_summary,
        final_k_summary,
    ) = compute_final_signatures(arms_by_disease, source, feature_rows, curve_all)
    final_features.to_csv(OUTDIR / "figure1_final_signature_features.tsv", sep="\t", index=False)
    final_auc.to_csv(OUTDIR / "figure1_final_signature_auc_by_cohort.tsv", sep="\t", index=False)
    final_predictions.to_csv(OUTDIR / "figure1_final_signature_predictions.tsv", sep="\t", index=False)
    final_disease_summary.to_csv(OUTDIR / "figure1_final_signature_disease_summary.tsv", sep="\t", index=False)
    final_k_summary.to_csv(OUTDIR / "figure1_final_signature_k_selection.tsv", sep="\t", index=False)

    exclusion_columns = ["disease", "study_code", "arm_id", "status", "message"]
    pd.DataFrame(exclusions, columns=exclusion_columns).to_csv(OUTDIR / "figure1_exclusions.tsv", sep="\t", index=False)

    if source.empty:
        raise RuntimeError("No LODO AUC rows were generated.")
    if not source["auc"].dropna().between(0, 1).all():
        raise RuntimeError("Generated AUC values outside [0, 1].")
    if final_auc.empty:
        raise RuntimeError("No final signature AUC rows were generated.")
    if not final_auc["auc"].dropna().between(0, 1).all():
        raise RuntimeError("Generated final signature AUC values outside [0, 1].")
    if not ((final_auc["auc_ci_low"] <= final_auc["auc"]) & (final_auc["auc"] <= final_auc["auc_ci_high"])).all():
        raise RuntimeError("Final signature AUC confidence intervals do not contain point estimates.")
    feature_counts = final_features.groupby("disease")["species"].nunique()
    effective_counts = final_auc.groupby("disease")["effective_k"].first()
    if not feature_counts.equals(effective_counts.reindex(feature_counts.index)):
        raise RuntimeError("Final signature feature counts do not match effective_k values.")
    print(f"Wrote {source.shape[0]} LODO rows to {OUTDIR / 'figure1_source_data.tsv'}", flush=True)


if __name__ == "__main__":
    run()
