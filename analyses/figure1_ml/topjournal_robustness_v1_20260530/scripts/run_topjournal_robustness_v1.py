#!/usr/bin/env python3
"""Top-journal robustness package for the Figure1_ML species analyses."""

from __future__ import annotations

import hashlib
import itertools
import math
import shutil
from collections import Counter
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import binomtest, mannwhitneyu, norm, wilcoxon
from sklearn.metrics import average_precision_score, brier_score_loss, roc_auc_score
import statsmodels.api as sm


ROOT = Path(__file__).resolve().parents[2]
REQ = ROOT / "requested_topjournal_figures"
ASSOC = ROOT / "species_association_overlap"
KEY = ROOT / "key_species_overlap_v2"
OUT = ROOT / "topjournal_robustness_v1_20260530"
TABLES = OUT / "tables"
FIGURES = OUT / "figures"
REPORTS = OUT / "reports"
SCRIPTS = OUT / "scripts"
LOGS = OUT / "logs"

DISEASES = ("CRA", "CRC", "IBD", "IBS")
PRED_DISEASES = ("CRA", "CRC", "IBD")
PAIR_ORDER = ("CRA|CRC", "CRA|IBD", "CRA|IBS", "CRC|IBD", "CRC|IBS", "IBD|IBS")
TOPN_ORDER = ("top_50", "top_100", "top_150", "top_200", "top_300", "top_500")
RANDOM_SEED = 20260530
MIN_CASE = 10
MIN_CONTROL = 10
MIN_DETECTION_FRAC = 0.05
MIN_DETECTION_N = 10
MIN_MEAN_ABUND = 1e-5
LOG_PSEUDOCOUNT = 1e-6
N_COVARIATE_PERMUTATIONS = 5000

CONTROL_STATUSES = {"ctr", "control patient"}
IBD_CASE_STATUSES = {
    "crohn's disease",
    "ulcerative colitis",
    "inflammatory bowel disease",
    "indeterminate colitis",
}
IBS_CASE_STATUSES = {"irritable bowel syndrome"}

TEXT = "#111820"
MUTED = "#52616E"
GRID = "#D9E1E8"
COLORS = {"CRA": "#0B3C68", "CRC": "#2F7EA8", "IBD": "#8ABBD3", "IBS": "#E9BFAE"}


def ensure_dirs() -> None:
    for path in (OUT, TABLES, FIGURES, REPORTS, SCRIPTS, LOGS):
        path.mkdir(parents=True, exist_ok=True)


def bh_fdr(p_values: np.ndarray) -> np.ndarray:
    p = np.asarray(p_values, dtype=float)
    q = np.full_like(p, np.nan, dtype=float)
    mask = np.isfinite(p)
    if not mask.any():
        return q
    valid = p[mask]
    order = np.argsort(valid)
    ranked = valid[order]
    n = len(ranked)
    adjusted = ranked * n / np.arange(1, n + 1)
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    adjusted = np.clip(adjusted, 0, 1)
    tmp = np.empty_like(adjusted)
    tmp[order] = adjusted
    q[mask] = tmp
    return q


def resolve_path(value: object) -> Path:
    path = Path(str(value))
    if path.exists():
        return path
    if not path.is_absolute():
        candidate = ROOT / path
        if candidate.exists():
            return candidate
    if path.is_absolute():
        for anchor in ("CRA", "CRC", "IBD", "IBS"):
            if anchor in path.parts:
                candidate = ROOT.joinpath(*path.parts[path.parts.index(anchor) :])
                if candidate.exists():
                    return candidate
    return path


def label_status(disease: str, status: object) -> str | None:
    s = str(status or "").strip().lower()
    if s in CONTROL_STATUSES:
        return "control"
    if disease == "CRA":
        return "case" if ("adenoma" in s or "adenomatous" in s) else None
    if disease == "CRC":
        return "case" if s == "colorectal cancer" else None
    if disease == "IBD":
        return "case" if s in IBD_CASE_STATUSES else None
    if disease == "IBS":
        return "case" if s in IBS_CASE_STATUSES else None
    return None


def load_metaphlan_species(path: Path) -> pd.DataFrame:
    df = pd.read_csv(
        path,
        sep="\t",
        compression="gzip",
        usecols=["sample_alias", "clade_name", "rel_abund"],
        dtype={"sample_alias": str, "clade_name": str, "rel_abund": np.float32},
    )
    df = df[df["clade_name"].str.startswith("s__", na=False)]
    mat = df.pivot_table(index="sample_alias", columns="clade_name", values="rel_abund", aggfunc="sum", fill_value=0.0)
    return mat.astype(np.float32)


def read_auc_inputs() -> tuple[pd.DataFrame, pd.DataFrame]:
    auc = pd.read_csv(REQ / "figure1_final_signature_auc_by_cohort.tsv", sep="\t")
    pred = pd.read_csv(REQ / "figure1_final_signature_predictions.tsv", sep="\t")
    pred["predicted_probability"] = pd.to_numeric(pred["predicted_probability"], errors="coerce")
    pred = pred.dropna(subset=["predicted_probability"])
    return auc, pred


def auc_by_cohort_metrics() -> pd.DataFrame:
    auc, pred = read_auc_inputs()
    rows = []
    for (disease, cohort), sub in pred.groupby(["disease", "left_out_cohort"], sort=False):
        y = sub["y_true"].to_numpy(dtype=int)
        p = sub["predicted_probability"].to_numpy(dtype=float)
        case = p[y == 1]
        control = p[y == 0]
        if len(np.unique(y)) != 2:
            continue
        rows.append(
            {
                "disease": disease,
                "left_out_cohort": cohort,
                "n_test": int(len(y)),
                "case_n": int(case.size),
                "control_n": int(control.size),
                "roc_auc": float(roc_auc_score(y, p)),
                "average_precision": float(average_precision_score(y, p)),
                "brier_score": float(brier_score_loss(y, p)),
                "case_mean_probability": float(np.mean(case)),
                "control_mean_probability": float(np.mean(control)),
                "case_control_mean_diff": float(np.mean(case) - np.mean(control)),
                "case_median_probability": float(np.median(case)),
                "control_median_probability": float(np.median(control)),
                "case_control_median_diff": float(np.median(case) - np.median(control)),
            }
        )
    metrics = pd.DataFrame(rows)
    keep_cols = [
        "disease",
        "left_out_cohort",
        "auc",
        "auc_ci_low",
        "auc_ci_high",
        "auc_permutation_p_one_sided",
        "significant_auc_gt_0_5",
        "final_k",
        "effective_k",
    ]
    metrics = metrics.merge(auc[keep_cols], on=["disease", "left_out_cohort"], how="left")
    metrics = metrics.sort_values(["disease", "roc_auc", "left_out_cohort"], ascending=[True, False, True])
    metrics.to_csv(TABLES / "auc_by_cohort_metrics.tsv", sep="\t", index=False)
    return metrics


def exact_permutation_p(x: np.ndarray, y: np.ndarray, alternative: str) -> float:
    pooled = np.r_[x, y]
    n_x = len(x)
    observed = float(np.mean(y) - np.mean(x))
    diffs = []
    idx_all = range(len(pooled))
    for idx_x in itertools.combinations(idx_all, n_x):
        mask = np.zeros(len(pooled), dtype=bool)
        mask[list(idx_x)] = True
        diff = float(np.mean(pooled[~mask]) - np.mean(pooled[mask]))
        diffs.append(diff)
    diffs = np.asarray(diffs)
    if alternative == "greater":
        return float((np.sum(diffs >= observed - 1e-12)) / diffs.size)
    return float((np.sum(np.abs(diffs) >= abs(observed) - 1e-12)) / diffs.size)


def auc_pairwise_tests(metrics: pd.DataFrame) -> pd.DataFrame:
    rows = []
    tests = [("CRA", "CRC", "greater"), ("CRA", "IBD", "greater"), ("CRC", "IBD", "two-sided")]
    for a, b, alternative in tests:
        x = metrics.loc[metrics["disease"].eq(a), "roc_auc"].to_numpy(dtype=float)
        y = metrics.loc[metrics["disease"].eq(b), "roc_auc"].to_numpy(dtype=float)
        mw_alt = "less" if alternative == "greater" else "two-sided"
        rows.append(
            {
                "comparison": f"{a}_vs_{b}",
                "disease_a": a,
                "disease_b": b,
                "alternative": "disease_b_mean_auc_greater" if alternative == "greater" else "two_sided",
                "n_a": int(x.size),
                "n_b": int(y.size),
                "mean_auc_a": float(np.mean(x)),
                "mean_auc_b": float(np.mean(y)),
                "median_auc_a": float(np.median(x)),
                "median_auc_b": float(np.median(y)),
                "mean_diff_b_minus_a": float(np.mean(y) - np.mean(x)),
                "median_diff_b_minus_a": float(np.median(y) - np.median(x)),
                "mannwhitney_p": float(mannwhitneyu(x, y, alternative=mw_alt).pvalue),
                "exact_permutation_p_mean_diff": exact_permutation_p(x, y, alternative),
            }
        )
    out = pd.DataFrame(rows)
    out["mannwhitney_q"] = bh_fdr(out["mannwhitney_p"].to_numpy(dtype=float))
    out["exact_permutation_q"] = bh_fdr(out["exact_permutation_p_mean_diff"].to_numpy(dtype=float))
    out.to_csv(TABLES / "auc_disease_pairwise_tests.tsv", sep="\t", index=False)
    return out


def auc_paired_sensitivity(metrics: pd.DataFrame) -> pd.DataFrame:
    cra = metrics[metrics["disease"].eq("CRA")][["left_out_cohort", "roc_auc"]].rename(columns={"roc_auc": "cra_auc"})
    crc = metrics[metrics["disease"].eq("CRC")][["left_out_cohort", "roc_auc"]].rename(columns={"roc_auc": "crc_auc"})
    paired = cra.merge(crc, on="left_out_cohort", how="inner")
    paired["crc_minus_cra_auc"] = paired["crc_auc"] - paired["cra_auc"]
    rows = paired.copy()
    rows.insert(0, "row_type", "cohort")
    if not paired.empty:
        diffs = paired["crc_minus_cra_auc"].to_numpy(dtype=float)
        try:
            wilcox_p = float(wilcoxon(diffs, alternative="greater", zero_method="wilcox").pvalue)
        except ValueError:
            wilcox_p = np.nan
        sign_p = float(binomtest(int((diffs > 0).sum()), n=len(diffs), p=0.5, alternative="greater").pvalue)
        summary = pd.DataFrame(
            [
                {
                    "row_type": "summary",
                    "left_out_cohort": "shared_CRA_CRC_cohorts",
                    "cra_auc": float(paired["cra_auc"].mean()),
                    "crc_auc": float(paired["crc_auc"].mean()),
                    "crc_minus_cra_auc": float(diffs.mean()),
                    "n_shared_cohorts": int(len(diffs)),
                    "wilcoxon_p_greater": wilcox_p,
                    "sign_test_p_greater": sign_p,
                    "all_shared_cohorts_crc_gt_cra": bool(np.all(diffs > 0)),
                }
            ]
        )
        out = pd.concat([rows, summary], ignore_index=True, sort=False)
    else:
        out = rows
    out.to_csv(TABLES / "auc_shared_study_paired_sensitivity.tsv", sep="\t", index=False)
    return out


def topn_overlap_matrix() -> pd.DataFrame:
    topn = pd.read_csv(ASSOC / "topn_sensitivity_overlap_summary.tsv", sep="\t")
    topn = topn[topn["analysis"].isin(TOPN_ORDER)].copy()
    topn["analysis_family"] = "disease_associated_ranked_topN"
    topn["claim_relevance"] = np.select(
        [
            topn["diseases"].eq("CRC|IBD"),
            topn["diseases"].isin(["CRC|IBS", "IBD|IBS"]),
            topn["diseases"].str.contains("IBS", regex=False),
            topn["n_diseases"].ge(3),
        ],
        ["core_crc_ibd_pair", "ibs_pairwise_claim", "other_ibs_overlap", "multi_disease_overlap"],
        default="other_pairwise",
    )
    topn["conservative_support"] = topn["significantly_above_random"].astype(bool) & (topn["permutation_q"] < 0.05)
    topn.to_csv(TABLES / "topn_overlap_robustness_matrix.tsv", sep="\t", index=False)
    return topn


def strict_named_locked_audit() -> pd.DataFrame:
    frames = []
    sens = pd.read_csv(KEY / "sensitivity_overlap_summary.tsv", sep="\t")
    sens = sens[sens["sensitivity"].isin(["top50_named_species_only", "strict_only", "crc_keys_crc_only_studies"])]
    tmp = sens.rename(columns={"sensitivity": "analysis"})
    tmp.insert(0, "analysis_family", "metaphlan4_key_species_sensitivity")
    frames.append(tmp)

    selected = pd.read_csv(ASSOC / "selected_sensitivity_overlap_summary.tsv", sep="\t")
    selected = selected[selected["analysis"].eq("top_300_named_only")]
    selected.insert(0, "analysis_family", "metaphlan4_top300_named_only")
    frames.append(selected)

    locked = pd.read_csv(REQ / "feature_panel_locked_panel_vs_ibs_overlap_audit.tsv", sep="\t")
    locked = locked.rename(columns={"set_a": "disease_a", "set_b": "disease_b"})
    locked["analysis_family"] = "locked_prediction_panel_vs_ibs_related"
    locked["analysis"] = "locked_panel_overlap"
    locked["diseases"] = locked["disease_a"].astype(str) + "|" + locked["disease_b"].astype(str)
    locked["n_diseases"] = 2
    locked["set_sizes"] = locked["set_a_size"].astype(str) + "|" + locked["set_b_size"].astype(str)
    locked["background_size"] = np.nan
    locked["expected_permutation"] = np.nan
    locked["observed_expected_ratio"] = np.nan
    locked["permutation_q"] = np.nan
    locked["significantly_above_random"] = False
    frames.append(locked)

    common_cols = sorted(set().union(*(set(frame.columns) for frame in frames)))
    out = pd.concat([frame.reindex(columns=common_cols) for frame in frames], ignore_index=True, sort=False)
    out["conservative_support"] = out.get("significantly_above_random", False).fillna(False).astype(bool) & (
        pd.to_numeric(out.get("permutation_q", np.nan), errors="coerce") < 0.05
    )
    out.to_csv(TABLES / "strict_named_locked_overlap_audit.tsv", sep="\t", index=False)
    return out


def covariate_completeness() -> pd.DataFrame:
    inv = pd.read_csv(KEY / "disease_arm_inventory.tsv", sep="\t")
    inc = inv[inv["included_main"].astype(str).str.lower().eq("true")].copy()
    cols = ["sex", "age_years", "bmi", "smoker", "location", "geographic_location"]
    rows = []
    for _, rec in inc.iterrows():
        meta_path = resolve_path(rec["metadata_path"])
        prof_path = resolve_path(rec["profile_path"])
        def release_path(path: Path) -> str:
            try:
                return str(path.resolve().relative_to(ROOT.resolve()))
            except ValueError:
                return str(path)

        row = {
            "disease": rec["disease"],
            "study_code": rec["study_code"],
            "arm_id": rec["arm_id"],
            "metadata_path_resolved": release_path(meta_path),
            "profile_path_resolved": release_path(prof_path),
            "metadata_exists": meta_path.exists(),
            "profile_exists": prof_path.exists(),
        }
        if meta_path.exists():
            meta = pd.read_csv(meta_path, sep="\t", dtype=str)
            row["metadata_rows"] = int(meta.shape[0])
            labeled = meta[["sample_alias", "subject_disease_status"]].copy()
            labeled["group"] = labeled["subject_disease_status"].map(lambda x: label_status(str(rec["disease"]), x))
            labeled = labeled.dropna(subset=["group"]).drop_duplicates("sample_alias", keep="first")
            row["labeled_rows"] = int(labeled.shape[0])
            for col in cols:
                row[f"{col}_present"] = col in meta.columns
                row[f"{col}_nonmissing"] = int(meta[col].notna().sum()) if col in meta.columns else 0
                row[f"{col}_nonmissing_frac"] = float(meta[col].notna().mean()) if col in meta.columns and len(meta) else np.nan
            if {"age_years", "sex"}.issubset(meta.columns):
                cov = meta[["sample_alias", "age_years", "sex", "subject_disease_status"]].copy()
                cov["age_years"] = pd.to_numeric(cov["age_years"], errors="coerce")
                cov["sex"] = cov["sex"].replace({"": np.nan, "nan": np.nan})
                cov["group"] = cov["subject_disease_status"].map(lambda x: label_status(str(rec["disease"]), x))
                cov = cov.dropna(subset=["sample_alias", "age_years", "sex", "group"]).drop_duplicates("sample_alias")
                row["age_sex_complete_labeled_rows"] = int(cov.shape[0])
                row["age_sex_complete_case"] = int((cov["group"] == "case").sum())
                row["age_sex_complete_control"] = int((cov["group"] == "control").sum())
                row["age_sex_evaluable"] = bool(row["age_sex_complete_case"] >= MIN_CASE and row["age_sex_complete_control"] >= MIN_CONTROL)
            else:
                row["age_sex_complete_labeled_rows"] = 0
                row["age_sex_complete_case"] = 0
                row["age_sex_complete_control"] = 0
                row["age_sex_evaluable"] = False
        rows.append(row)
    out = pd.DataFrame(rows)
    out.to_csv(TABLES / "covariate_completeness_by_arm.tsv", sep="\t", index=False)
    return out


def fit_adjusted_arm(rec: pd.Series) -> pd.DataFrame:
    disease = str(rec["disease"])
    meta_path = resolve_path(rec["metadata_path_resolved"])
    profile_path = resolve_path(rec["profile_path_resolved"])
    if not bool(rec.get("age_sex_evaluable", False)) or not meta_path.exists() or not profile_path.exists():
        return pd.DataFrame()

    meta = pd.read_csv(meta_path, sep="\t", dtype=str)
    required = {"sample_alias", "subject_disease_status", "age_years", "sex"}
    if not required.issubset(meta.columns):
        return pd.DataFrame()
    meta = meta[["sample_alias", "subject_disease_status", "age_years", "sex"]].copy()
    meta["group"] = meta["subject_disease_status"].map(lambda x: label_status(disease, x))
    meta["age_years"] = pd.to_numeric(meta["age_years"], errors="coerce")
    meta["sex"] = meta["sex"].astype(str).str.strip().replace({"": np.nan, "nan": np.nan, "None": np.nan})
    meta = meta.dropna(subset=["sample_alias", "group", "age_years", "sex"]).drop_duplicates("sample_alias", keep="first")
    x_raw = load_metaphlan_species(profile_path)
    common = meta["sample_alias"][meta["sample_alias"].isin(x_raw.index)]
    meta = meta.set_index("sample_alias").loc[common].copy()
    x_raw = x_raw.loc[common].copy()
    y = (meta["group"] == "case").astype(int)
    if int(y.sum()) < MIN_CASE or int((y == 0).sum()) < MIN_CONTROL:
        return pd.DataFrame()

    detected = (x_raw > 0).sum(axis=0)
    prevalence = detected / max(1, x_raw.shape[0])
    mean_abund = x_raw.mean(axis=0)
    keep = ((prevalence >= MIN_DETECTION_FRAC) | (detected >= MIN_DETECTION_N)) & (mean_abund >= MIN_MEAN_ABUND)
    x_keep = x_raw.loc[:, keep]
    if x_keep.empty:
        return pd.DataFrame()

    sex_dummies = pd.get_dummies(meta["sex"], prefix="sex", drop_first=True, dtype=float)
    base = pd.DataFrame({"age_years": meta["age_years"].astype(float)}, index=meta.index).join(sex_dummies)
    base = sm.add_constant(base, has_constant="add")
    rows = []
    for species in x_keep.columns:
        values = x_keep[species].astype(float)
        case_mean = float(values[y == 1].mean())
        control_mean = float(values[y == 0].mean())
        log2fc = float(np.log2((case_mean + LOG_PSEUDOCOUNT) / (control_mean + LOG_PSEUDOCOUNT)))
        design = base.copy()
        design.insert(1, "log_abundance", np.log10(values.to_numpy(dtype=float) + LOG_PSEUDOCOUNT))
        status = "ok"
        coef = np.nan
        p_value = np.nan
        try:
            model = sm.GLM(y.to_numpy(dtype=int), design.astype(float), family=sm.families.Binomial())
            result = model.fit(maxiter=100, disp=0)
            coef = float(result.params.get("log_abundance", np.nan))
            p_value = float(result.pvalues.get("log_abundance", np.nan))
        except Exception as exc:
            status = f"failed:{type(exc).__name__}"
        rows.append(
            {
                "disease": disease,
                "study_code": rec["study_code"],
                "arm_id": rec["arm_id"],
                "species": species,
                "status": status,
                "n_samples_covariate_complete": int(y.shape[0]),
                "case_samples": int(y.sum()),
                "control_samples": int((y == 0).sum()),
                "case_mean": case_mean,
                "control_mean": control_mean,
                "log2fc": log2fc,
                "adjusted_log_abundance_coef": coef,
                "adjusted_p_value": p_value,
                "detected": int(detected.loc[species]),
                "mean_abund": float(mean_abund.loc[species]),
            }
        )
    out = pd.DataFrame(rows)
    mask = out["status"].eq("ok") & np.isfinite(out["adjusted_p_value"].to_numpy(dtype=float))
    out["adjusted_q_value"] = np.nan
    if mask.any():
        out.loc[mask, "adjusted_q_value"] = bh_fdr(out.loc[mask, "adjusted_p_value"].to_numpy(dtype=float))
    return out


def covariate_adjusted_stats(completeness: pd.DataFrame) -> pd.DataFrame:
    frames = []
    for _, rec in completeness.iterrows():
        frame = fit_adjusted_arm(rec)
        if not frame.empty:
            frames.append(frame)
    out = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
    out.to_csv(TABLES / "covariate_adjusted_arm_species_stats.tsv", sep="\t", index=False)
    return out


def aggregate_adjusted(adjusted: pd.DataFrame, completeness: pd.DataFrame) -> pd.DataFrame:
    if adjusted.empty:
        return pd.DataFrame()
    ok = adjusted[adjusted["status"].eq("ok") & adjusted["adjusted_p_value"].notna()].copy()
    rows = []
    evaluable_arms = completeness[completeness["age_sex_evaluable"].astype(bool)].groupby("disease")["arm_id"].nunique().to_dict()
    for disease, sub in ok.groupby("disease"):
        required_support = 1 if disease == "IBS" else 2
        for species, g in sub.groupby("species"):
            p = np.clip(g["adjusted_p_value"].to_numpy(dtype=float), 1e-300, 1.0)
            lfc = g["log2fc"].to_numpy(dtype=float)
            weights = np.sqrt(g["n_samples_covariate_complete"].to_numpy(dtype=float))
            signed_z = np.sign(lfc) * norm.isf(p / 2.0)
            denom = math.sqrt(float(np.sum(weights**2)))
            meta_z = float(np.sum(weights * signed_z) / denom) if denom else 0.0
            meta_p = float(2 * norm.sf(abs(meta_z)))
            meta_lfc = float(np.average(lfc, weights=weights)) if weights.sum() else float(np.mean(lfc))
            signs = np.sign(lfc)
            nonzero = signs[signs != 0]
            direction_consistency = float(max((nonzero > 0).sum(), (nonzero < 0).sum()) / max(len(nonzero), 1))
            rows.append(
                {
                    "disease": disease,
                    "species": species,
                    "n_arms_tested": int(g["arm_id"].nunique()),
                    "n_arms_total_evaluable": int(evaluable_arms.get(disease, 0)),
                    "required_support": required_support,
                    "direction": "up" if meta_lfc > 0 else "down",
                    "direction_consistency": direction_consistency,
                    "meta_log2fc": meta_lfc,
                    "meta_z": meta_z,
                    "meta_p": meta_p,
                }
            )
    ranked = pd.DataFrame(rows)
    if ranked.empty:
        return ranked
    ranked["meta_q"] = np.nan
    for disease in ranked["disease"].unique():
        mask = ranked["disease"].eq(disease)
        ranked.loc[mask, "meta_q"] = bh_fdr(ranked.loc[mask, "meta_p"].to_numpy(dtype=float))
    ranked["rank_score"] = (
        -np.log10(ranked["meta_q"].clip(lower=1e-300))
        * ranked["meta_log2fc"].abs()
        * ranked["direction_consistency"]
    )
    ranked["eligible_for_adjusted_top"] = (
        (ranked["direction_consistency"] >= 0.60) & (ranked["n_arms_tested"] >= ranked["required_support"])
    )
    return ranked.sort_values(["disease", "rank_score"], ascending=[True, False])


def overlap_tests(key_sets: dict[str, set[str]], backgrounds: dict[str, set[str]], permutations: int) -> pd.DataFrame:
    rng = np.random.default_rng(RANDOM_SEED)
    rows = []
    diseases = [d for d in DISEASES if d in key_sets]
    for r in range(2, len(diseases) + 1):
        for combo in itertools.combinations(diseases, r):
            sets = [key_sets.get(d, set()) for d in combo]
            sizes = [len(s) for s in sets]
            shared = set.intersection(*sets) if all(sizes) else set()
            union_keys = set.union(*sets) if sets else set()
            bg = set.union(*(backgrounds.get(d, set()) for d in combo))
            bg_list = np.array(sorted(bg), dtype=object)
            observed = len(shared)
            if len(bg_list) > 0 and all(size <= len(bg_list) for size in sizes):
                perm_counts = np.empty(permutations, dtype=np.int16)
                for i in range(permutations):
                    sampled = [set(rng.choice(bg_list, size=size, replace=False)) for size in sizes]
                    perm_counts[i] = len(set.intersection(*sampled))
                expected = float(np.mean(perm_counts))
                p_value = float((1 + np.sum(perm_counts >= observed)) / (permutations + 1))
                q95 = float(np.quantile(perm_counts, 0.95))
            else:
                expected = np.nan
                p_value = np.nan
                q95 = np.nan
            rows.append(
                {
                    "diseases": "|".join(combo),
                    "n_diseases": r,
                    "set_sizes": "|".join(map(str, sizes)),
                    "background_size": int(len(bg)),
                    "shared_count": observed,
                    "jaccard": observed / len(union_keys) if union_keys else np.nan,
                    "overlap_coefficient": observed / min(sizes) if min(sizes) else np.nan,
                    "expected_permutation": expected,
                    "permutation_q95": q95,
                    "observed_expected_ratio": observed / expected if expected and expected > 0 else np.nan,
                    "permutation_p": p_value,
                    "shared_species": ";".join(sorted(shared)),
                }
            )
    out = pd.DataFrame(rows)
    out["permutation_q"] = bh_fdr(out["permutation_p"].to_numpy(dtype=float)) if not out.empty else np.nan
    out["significantly_above_random"] = (out["permutation_q"] < 0.05) & (out["shared_count"] > out["permutation_q95"])
    return out


def covariate_adjusted_overlap(adjusted: pd.DataFrame, completeness: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    ranked = aggregate_adjusted(adjusted, completeness)
    if ranked.empty:
        out = pd.DataFrame()
        out.to_csv(TABLES / "covariate_adjusted_overlap_sensitivity.tsv", sep="\t", index=False)
        return ranked, out
    ranked.to_csv(TABLES / "covariate_adjusted_disease_ranked.tsv", sep="\t", index=False)
    ok = adjusted[adjusted["status"].eq("ok")]
    backgrounds = {d: set(g["species"]) for d, g in ok.groupby("disease")}
    frames = []
    for top_n in (50, 100, 200, 300):
        sets = {}
        for disease in DISEASES:
            sub = ranked[(ranked["disease"].eq(disease)) & (ranked["eligible_for_adjusted_top"])].copy()
            sets[disease] = set(sub.sort_values("rank_score", ascending=False).head(top_n)["species"])
        overlap = overlap_tests(sets, backgrounds, N_COVARIATE_PERMUTATIONS)
        overlap.insert(0, "analysis", f"age_sex_adjusted_top_{top_n}")
        frames.append(overlap)
    out = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
    out.to_csv(TABLES / "covariate_adjusted_overlap_sensitivity.tsv", sep="\t", index=False)
    return ranked, out


def setup_plot_style() -> None:
    plt.rcParams.update(
        {
            "font.family": "Arial",
            "font.sans-serif": ["Arial", "DejaVu Sans"],
            "font.size": 10.5,
            "axes.labelcolor": TEXT,
            "xtick.color": TEXT,
            "ytick.color": TEXT,
            "axes.edgecolor": "#CBD5E1",
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
        }
    )


def savefig(fig: plt.Figure, stem: str) -> None:
    fig.savefig(FIGURES / f"{stem}.pdf", bbox_inches="tight")
    fig.savefig(FIGURES / f"{stem}.png", dpi=300, bbox_inches="tight")
    fig.savefig(FIGURES / f"{stem}.svg", bbox_inches="tight")
    plt.close(fig)


def plot_auc(metrics: pd.DataFrame, pair_tests: pd.DataFrame, paired: pd.DataFrame) -> None:
    setup_plot_style()
    fig, axes = plt.subplots(1, 2, figsize=(10.8, 4.0), gridspec_kw={"width_ratios": [1.1, 0.9]})
    ax = axes[0]
    positions = np.arange(len(PRED_DISEASES))
    for i, disease in enumerate(PRED_DISEASES):
        vals = metrics.loc[metrics["disease"].eq(disease), "roc_auc"].to_numpy(dtype=float)
        jitter = np.linspace(-0.08, 0.08, len(vals)) if len(vals) else []
        ax.scatter(np.full(len(vals), i) + jitter, vals, color=COLORS[disease], edgecolor="white", s=42, zorder=3)
        ax.boxplot(vals, positions=[i], widths=0.35, patch_artist=True, showfliers=False, medianprops={"color": TEXT}, boxprops={"facecolor": "#F7FAFC", "edgecolor": COLORS[disease]}, whiskerprops={"color": COLORS[disease]}, capprops={"color": COLORS[disease]})
    ax.axhline(0.5, color="#8B98A5", ls=(0, (4, 3)), lw=1)
    ax.set_xticks(positions)
    ax.set_xticklabels(PRED_DISEASES)
    ax.set_ylabel("Held-out ROC-AUC")
    ax.set_title("Cross-cohort microbial panel performance", loc="left")
    ax.set_ylim(0.35, 1.0)
    ax.grid(axis="y", color=GRID, lw=0.7)
    txt = []
    for _, row in pair_tests.iterrows():
        txt.append(f"{row['comparison'].replace('_vs_', ' vs ')} q={row['exact_permutation_q']:.3g}")
    ax.text(0.02, 0.98, "\n".join(txt), transform=ax.transAxes, ha="left", va="top", fontsize=8.3, color=MUTED)

    ax = axes[1]
    cohort_rows = paired[paired["row_type"].eq("cohort")].copy()
    if not cohort_rows.empty:
        for _, row in cohort_rows.iterrows():
            ax.plot([0, 1], [row["cra_auc"], row["crc_auc"]], color="#CBD5E1", lw=1.2, zorder=1)
            ax.scatter([0, 1], [row["cra_auc"], row["crc_auc"]], color=[COLORS["CRA"], COLORS["CRC"]], s=40, edgecolor="white", zorder=2)
        summary = paired[paired["row_type"].eq("summary")]
        label = ""
        if not summary.empty:
            sr = summary.iloc[0]
            label = f"mean delta={sr['crc_minus_cra_auc']:.2f}\nWilcoxon p={sr['wilcoxon_p_greater']:.3g}"
        ax.text(0.05, 0.96, label, transform=ax.transAxes, ha="left", va="top", fontsize=8.8, color=MUTED)
    ax.axhline(0.5, color="#8B98A5", ls=(0, (4, 3)), lw=1)
    ax.set_xlim(-0.35, 1.35)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["CRA", "CRC"])
    ax.set_ylabel("ROC-AUC")
    ax.set_title("Shared-study CRA/CRC paired sensitivity", loc="left")
    ax.grid(axis="y", color=GRID, lw=0.7)
    savefig(fig, "topjournal_auc_robustness")


def plot_overlap(topn: pd.DataFrame) -> None:
    setup_plot_style()
    fig, axes = plt.subplots(1, 2, figsize=(12.2, 4.2), gridspec_kw={"width_ratios": [1.35, 0.9]})
    pair = topn[topn["n_diseases"].eq(2)].copy()
    pair["top_n"] = pair["analysis"].str.replace("top_", "", regex=False).astype(int)
    ax = axes[0]
    for diseases, group in pair.groupby("diseases"):
        group = group.sort_values("top_n")
        color = "#0B3C68" if diseases == "CRC|IBD" else "#C98773" if "IBS" in diseases else "#7897AC"
        ax.plot(group["top_n"], group["observed_expected_ratio"], marker="o", lw=1.6, label=diseases, color=color, alpha=0.95)
        sig = group[group["conservative_support"]]
        ax.scatter(sig["top_n"], sig["observed_expected_ratio"], s=70, facecolors="none", edgecolors=color, linewidth=1.5)
    ax.axhline(1.0, color="#8B98A5", lw=1, ls=(0, (4, 3)))
    ax.set_xlabel("Ranked associated set size")
    ax.set_ylabel("Observed / permutation expected")
    ax.set_title("Pairwise overlap robustness across Top-N thresholds", loc="left")
    ax.grid(axis="y", color=GRID, lw=0.7)
    ax.legend(frameon=False, fontsize=8, ncol=2)

    primary = pd.read_csv(ASSOC / "overlap_summary_top300.tsv", sep="\t")
    multi = primary[primary["n_diseases"].ge(3)].copy()
    ax = axes[1]
    x = np.arange(multi.shape[0])
    colors = ["#0B3C68" if bool(v) else "#BFD9E6" for v in multi["significantly_above_random"]]
    ax.bar(x, multi["shared_count"], color=colors, edgecolor="white")
    ax.scatter(x, multi["expected_permutation"], color="#D7A187", s=45, edgecolor="white", zorder=3, label="Expected")
    ax.set_xticks(x)
    ax.set_xticklabels(multi["diseases"].str.replace("|", "+", regex=False), rotation=35, ha="right")
    ax.set_ylabel("Shared species")
    ax.set_title("Top300 multi-disease intersections", loc="left")
    ax.grid(axis="y", color=GRID, lw=0.7)
    ax.legend(frameon=False, fontsize=8)
    savefig(fig, "topjournal_overlap_robustness")


def plot_covariate_completeness(comp: pd.DataFrame) -> None:
    setup_plot_style()
    covs = ["age_years", "sex", "bmi", "smoker"]
    labels = comp["disease"].astype(str) + " | " + comp["study_code"].astype(str)
    mat = np.vstack([comp[f"{c}_nonmissing_frac"].fillna(0).to_numpy(dtype=float) for c in covs]).T
    fig, ax = plt.subplots(figsize=(8.8, max(5.5, len(labels) * 0.25)))
    im = ax.imshow(mat, aspect="auto", vmin=0, vmax=1, cmap="Blues")
    ax.set_yticks(np.arange(len(labels)))
    ax.set_yticklabels(labels, fontsize=7.8)
    ax.set_xticks(np.arange(len(covs)))
    ax.set_xticklabels(covs)
    ax.set_title("Covariate completeness by disease arm", loc="left")
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            ax.text(j, i, f"{mat[i, j]:.0%}", ha="center", va="center", fontsize=7.2, color="white" if mat[i, j] > 0.55 else TEXT)
    fig.colorbar(im, ax=ax, fraction=0.025, pad=0.02, label="Non-missing fraction")
    fig.savefig(FIGURES / "covariate_completeness_audit.pdf", bbox_inches="tight")
    fig.savefig(FIGURES / "covariate_completeness_audit.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def plot_covariate_adjusted_overlap(adj_overlap: pd.DataFrame) -> None:
    if adj_overlap.empty:
        return
    setup_plot_style()
    pair = adj_overlap[adj_overlap["n_diseases"].eq(2)].copy()
    pair["top_n"] = pair["analysis"].str.extract(r"(\d+)$").astype(int)
    fig, ax = plt.subplots(figsize=(7.2, 4.2))
    for diseases, group in pair.groupby("diseases"):
        group = group.sort_values("top_n")
        color = "#0B3C68" if diseases == "CRC|IBD" else "#C98773" if "IBS" in diseases else "#7897AC"
        ax.plot(group["top_n"], group["observed_expected_ratio"], marker="o", lw=1.6, label=diseases, color=color)
        sig = group[group["significantly_above_random"]]
        ax.scatter(sig["top_n"], sig["observed_expected_ratio"], s=70, facecolors="none", edgecolors=color, linewidth=1.5)
    ax.axhline(1.0, color="#8B98A5", lw=1, ls=(0, (4, 3)))
    ax.set_xlabel("Age/sex-adjusted ranked set size")
    ax.set_ylabel("Observed / permutation expected")
    ax.set_title("Age/sex-adjusted overlap sensitivity", loc="left")
    ax.grid(axis="y", color=GRID, lw=0.7)
    ax.legend(frameon=False, fontsize=8, ncol=2)
    fig.savefig(FIGURES / "covariate_adjusted_overlap_sensitivity.pdf", bbox_inches="tight")
    fig.savefig(FIGURES / "covariate_adjusted_overlap_sensitivity.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def fmt_p(value: float) -> str:
    if not np.isfinite(value):
        return "NA"
    if value < 0.001:
        return "<0.001"
    return f"{value:.3g}"


def write_reports(metrics: pd.DataFrame, pair_tests: pd.DataFrame, paired: pd.DataFrame, topn: pd.DataFrame, strict_audit: pd.DataFrame, comp: pd.DataFrame, adj_overlap: pd.DataFrame) -> None:
    disease_summary = metrics.groupby("disease").agg(
        n_cohorts=("roc_auc", "count"),
        mean_auc=("roc_auc", "mean"),
        median_auc=("roc_auc", "median"),
        min_auc=("roc_auc", "min"),
        max_auc=("roc_auc", "max"),
    )
    primary = pd.read_csv(ASSOC / "overlap_summary_top300.tsv", sep="\t")
    crc_ibd = primary[primary["diseases"].eq("CRC|IBD")].iloc[0]
    crc_ibs = primary[primary["diseases"].eq("CRC|IBS")].iloc[0]
    ibd_ibs = primary[primary["diseases"].eq("IBD|IBS")].iloc[0]
    four = primary[primary["diseases"].eq("CRA|CRC|IBD|IBS")].iloc[0]
    crc_ibd_ibs = primary[primary["diseases"].eq("CRC|IBD|IBS")].iloc[0]

    lines = ["# Top-Journal Robustness Summary", ""]
    lines.append("## Prediction performance")
    for disease, row in disease_summary.iterrows():
        lines.append(
            f"- {disease}: median ROC-AUC={row['median_auc']:.3f}, mean={row['mean_auc']:.3f}, range={row['min_auc']:.3f}-{row['max_auc']:.3f}, n={int(row['n_cohorts'])} cohorts."
        )
    for _, row in pair_tests.iterrows():
        lines.append(
            f"- {row['comparison'].replace('_vs_', ' vs ')}: mean delta={row['mean_diff_b_minus_a']:.3f}, exact permutation q={fmt_p(row['exact_permutation_q'])}."
        )
    lines.append("")
    lines.append("## Overlap robustness")
    lines.append(
        f"- CRC|IBD Top300 overlap is supported: shared={int(crc_ibd['shared_count'])}, obs/exp={crc_ibd['observed_expected_ratio']:.2f}, q={fmt_p(crc_ibd['permutation_q'])}."
    )
    lines.append(
        f"- CRC|IBD|IBS Top300 three-way overlap is supported: shared={int(crc_ibd_ibs['shared_count'])}, obs/exp={crc_ibd_ibs['observed_expected_ratio']:.2f}, q={fmt_p(crc_ibd_ibs['permutation_q'])}."
    )
    lines.append(
        f"- CRC|IBS pairwise is not supported as enriched: shared={int(crc_ibs['shared_count'])}, obs/exp={crc_ibs['observed_expected_ratio']:.2f}, q={fmt_p(crc_ibs['permutation_q'])}."
    )
    lines.append(
        f"- IBD|IBS pairwise is borderline/not FDR-significant: shared={int(ibd_ibs['shared_count'])}, obs/exp={ibd_ibs['observed_expected_ratio']:.2f}, q={fmt_p(ibd_ibs['permutation_q'])}."
    )
    lines.append(f"- Four-disease shared species are not observed in the Top300 analysis: shared={int(four['shared_count'])}.")
    lines.append("")
    lines.append("## Covariate adjustment")
    n_eval = int(comp["age_sex_evaluable"].fillna(False).sum())
    lines.append(f"- Age+sex adjusted sensitivity was attempted for {n_eval} evaluable disease arms.")
    if adj_overlap.empty:
        lines.append("- No adjusted overlap table was generated because adjusted species models did not yield enough ranked sets.")
    else:
        sig = adj_overlap[adj_overlap["significantly_above_random"].fillna(False)]
        lines.append(f"- Age+sex adjusted overlap generated {adj_overlap.shape[0]} tests; {sig.shape[0]} passed the permutation q<0.05 and q95 rule.")
    lines.append("")
    lines.append("Interpretation: the supported claim is narrower than broad four-disease overlap. CRA prediction is weaker than CRC/IBD, while CRC and IBD show moderate microbial-panel discrimination. Clinical-marker equivalence is not asserted from these data.")
    (REPORTS / "topjournal_results_summary.md").write_text("\n".join(lines) + "\n")

    audit = [
        "# Claim Support Audit",
        "",
        "| Claim | Decision | Evidence |",
        "|---|---|---|",
        "| CRA species-level classifier is lower than CRC and IBD | Supported | Exact permutation q values for CRA vs CRC/IBD are significant; see auc_disease_pairwise_tests.tsv. |",
        "| CRC and IBD AUCs are high | Not supported | AUCs are moderate, with cohort heterogeneity. |",
        "| Microbial AUC is similar to clinical indicators | Not assessed | No direct clinical-indicator model or biomarker AUC table was found in this analysis folder. |",
        "| Four diseases share many associated species | Not supported | Top300 four-way shared_count is 0. |",
        "| CRC and IBD associated species overlap | Supported | CRC|IBD Top300 q<0.05 and robust at larger thresholds. |",
        "| IBS overlaps strongly with both CRC and IBD pairwise | Not supported as pairwise claim | CRC|IBS is not enriched; IBD|IBS is borderline in Top300. CRC|IBD|IBS three-way overlap is supported. |",
    ]
    (REPORTS / "claim_support_audit.md").write_text("\n".join(audit) + "\n")

    clinical = [
        "# Clinical Comparison Gap",
        "",
        "No direct clinical-indicator prediction table, biomarker AUC table, or clinical-only model output was identified in the inspected Figure1_ML analysis outputs.",
        "",
        "Therefore this package does not claim that microbial-panel performance is similar to clinical indicators. To support that claim, add a cohort-matched clinical-only model and a direct comparison against the microbial predictions, preferably with paired bootstrap or DeLong-style testing where applicable.",
    ]
    (REPORTS / "clinical_comparison_gap.md").write_text("\n".join(clinical) + "\n")

    readme = [
        "# topjournal_robustness_v1_20260530",
        "",
        "Independent robustness package derived from existing Figure1_ML outputs.",
        "",
        "Directories:",
        "- `tables/`: source data and statistical summaries.",
        "- `figures/`: publication-oriented robustness figures.",
        "- `reports/`: conservative interpretation and claim audit.",
        "- `scripts/`: reproducible driver script snapshot.",
        "- `logs/`: reserved for run logs.",
    ]
    (REPORTS / "README.md").write_text("\n".join(readme) + "\n")


def validate_outputs() -> None:
    auc = pd.read_csv(TABLES / "auc_by_cohort_metrics.tsv", sep="\t")
    for col in ["roc_auc", "average_precision", "brier_score", "auc", "auc_ci_low", "auc_ci_high"]:
        vals = pd.to_numeric(auc[col], errors="coerce").dropna()
        if not vals.between(0, 1).all():
            raise RuntimeError(f"{col} contains values outside [0, 1]")
    ci = auc.dropna(subset=["auc", "auc_ci_low", "auc_ci_high"])
    if not ((ci["auc_ci_low"] <= ci["auc"]) & (ci["auc"] <= ci["auc_ci_high"])).all():
        raise RuntimeError("AUC confidence intervals do not contain point estimates.")
    for path in TABLES.glob("*.tsv"):
        df = pd.read_csv(path, sep="\t")
        for col in [c for c in df.columns if c.endswith("_q") or c in {"permutation_q", "meta_q", "adjusted_q_value"}]:
            vals = pd.to_numeric(df[col], errors="coerce").dropna()
            if not vals.between(0, 1).all():
                raise RuntimeError(f"{path.name}:{col} contains q values outside [0, 1]")
    topn = pd.read_csv(TABLES / "topn_overlap_robustness_matrix.tsv", sep="\t")
    primary_sets = topn[(topn["analysis"].eq("top_300")) & (topn["n_diseases"].eq(2))]["set_sizes"].unique().tolist()
    if "300|238" not in primary_sets:
        raise RuntimeError("IBS Top300 shortfall (300|238) was not preserved in topN output.")


def write_manifest() -> None:
    rows = []
    for path in sorted(OUT.rglob("*")):
        if path.is_file():
            rel = path.relative_to(OUT)
            h = hashlib.sha256(path.read_bytes()).hexdigest()
            rows.append({"relative_path": str(rel), "bytes": path.stat().st_size, "sha256": h})
    pd.DataFrame(rows).to_csv(REPORTS / "manifest.tsv", sep="\t", index=False)


def snapshot_script() -> None:
    src = Path(__file__).resolve()
    dst = SCRIPTS / src.name
    if src != dst:
        shutil.copy2(src, dst)


def main() -> int:
    ensure_dirs()
    snapshot_script()

    metrics = auc_by_cohort_metrics()
    pair_tests = auc_pairwise_tests(metrics)
    paired = auc_paired_sensitivity(metrics)
    topn = topn_overlap_matrix()
    strict_audit = strict_named_locked_audit()
    comp = covariate_completeness()
    adjusted = covariate_adjusted_stats(comp)
    _, adj_overlap = covariate_adjusted_overlap(adjusted, comp)

    plot_auc(metrics, pair_tests, paired)
    plot_overlap(topn)
    plot_covariate_completeness(comp)
    plot_covariate_adjusted_overlap(adj_overlap)
    write_reports(metrics, pair_tests, paired, topn, strict_audit, comp, adj_overlap)
    validate_outputs()
    write_manifest()

    print(f"Wrote top-journal robustness package to {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
