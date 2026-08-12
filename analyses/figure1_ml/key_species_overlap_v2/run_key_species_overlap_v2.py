#!/usr/bin/env python3
"""V2 key species overlap analysis using standardized Metalog disease folders."""

from __future__ import annotations

import argparse
import csv
import itertools
import math
import os
from collections import Counter
from dataclasses import dataclass, replace
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", str(Path(__file__).resolve().parent / ".mplconfig"))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import hypergeom, mannwhitneyu, norm
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold


ROOT = Path(__file__).resolve().parents[1]
OUTDIR = ROOT / "key_species_overlap_v2"
FIGDIR = OUTDIR / "figures"

DISEASES = ("CRA", "CRC", "IBD", "IBS")
CONTROL_STATUSES = {"ctr", "control patient"}
IBD_CASE_STATUSES = {
    "crohn's disease",
    "ulcerative colitis",
    "inflammatory bowel disease",
    "indeterminate colitis",
}
IBS_CASE_STATUSES = {"irritable bowel syndrome"}

MIN_CASE = 10
MIN_CONTROL = 10
MIN_DETECTION_FRAC = 0.05
MIN_DETECTION_N = 10
MIN_MEAN_ABUND = 1e-5
LOG_PSEUDOCOUNT = 1e-6
DIFF_Q_CUTOFF = 0.10
LOG2FC_CUTOFF = 0.50
DIRECTION_CONSISTENCY_CUTOFF = 0.60
TOP_N_MAIN = 50
RANDOM_SEED = 20260524

CRC_ONLY_STUDIES = {
    "Liu_2022_CRC_China",
    "Vogtmann_2016_CRC_USA",
    "Wirbel_2019_CRC_Germany",
    "Yu_2017_CRC_China",
}


@dataclass(frozen=True)
class ArmConfig:
    disease: str
    source_disease: str
    study_code: str
    arm_id: str
    metadata_path: Path | None
    profile_path: Path | None
    comparison: str
    skip_reason: str = ""


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


def find_profile(study_dir: Path, prefix: str = "metaphlan4_species") -> Path | None:
    latest = sorted(study_dir.glob(f"{prefix}_*_latest.tsv.gz"))
    if latest:
        return latest[0]
    legacy = sorted(study_dir.glob(f"{prefix}*.tsv.gz"))
    return legacy[0] if legacy else None


def profile_is_readable(path: Path | None) -> bool:
    if path is None or not path.exists() or path.stat().st_size == 0:
        return False
    try:
        pd.read_csv(
            path,
            sep="\t",
            compression="gzip",
            usecols=["sample_alias", "clade_name", "rel_abund"],
            nrows=5,
        )
        return True
    except Exception:
        return False


def metadata_counts(path: Path | None) -> str:
    if path is None or not path.exists():
        return ""
    counts: Counter[str] = Counter()
    try:
        with path.open("r", newline="") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                label = (row.get("subject_disease_status") or "").strip()
                counts[label or "<blank>"] += 1
    except Exception:
        return ""
    return "; ".join(f"{k}:{counts[k]}" for k in sorted(counts))


def make_config(disease: str, source_disease: str, study_dir: Path, comparison: str) -> ArmConfig:
    study_code = study_dir.name
    return ArmConfig(
        disease=disease,
        source_disease=source_disease,
        study_code=study_code,
        arm_id=f"{disease}__{source_disease}__{study_code}",
        metadata_path=study_dir / "metadata_all_wide.tsv",
        profile_path=find_profile(study_dir),
        comparison=comparison,
    )


def preferred_config(configs: list[ArmConfig]) -> ArmConfig:
    def score(cfg: ArmConfig) -> tuple[int, int, int]:
        return (
            1 if cfg.source_disease == "CRC" else 0,
            1 if profile_is_readable(cfg.profile_path) else 0,
            1 if cfg.metadata_path and cfg.metadata_path.exists() else 0,
        )

    return sorted(configs, key=score, reverse=True)[0]


def build_arm_configs() -> list[ArmConfig]:
    candidates: list[ArmConfig] = []

    for study_dir in sorted((ROOT / "CRA").glob("*")):
        if study_dir.is_dir() and (study_dir / "metadata_all_wide.tsv").exists():
            candidates.append(
                make_config(
                    "CRA",
                    "CRA",
                    study_dir,
                    "CRA: adenoma/adenomatous status vs CTR/control patient",
                )
            )

    for study_dir in sorted((ROOT / "CRC").glob("*")):
        if study_dir.is_dir() and (study_dir / "metadata_all_wide.tsv").exists():
            candidates.append(
                make_config(
                    "CRA",
                    "CRC",
                    study_dir,
                    "CRA arm extracted from CRC study: adenoma status vs CTR/control patient",
                )
            )
            candidates.append(
                make_config(
                    "CRC",
                    "CRC",
                    study_dir,
                    "CRC: colorectal cancer vs CTR/control patient; adenoma samples excluded",
                )
            )

    for source in ("IBD", "IBS"):
        for study_dir in sorted((ROOT / source).glob("*")):
            if not study_dir.is_dir() or not (study_dir / "metadata_all_wide.tsv").exists():
                continue
            if source == "IBS" and study_dir.name == "Mars_2020":
                continue
            comparison = (
                "IBD: Crohn/UC/IBD/indeterminate colitis vs CTR/control patient"
                if source == "IBD"
                else "IBS: irritable bowel syndrome vs CTR/control patient"
            )
            candidates.append(make_config(source, source, study_dir, comparison))

    grouped: dict[tuple[str, str], list[ArmConfig]] = {}
    for cfg in candidates:
        grouped.setdefault((cfg.disease, cfg.study_code), []).append(cfg)

    selected: list[ArmConfig] = []
    for _, group in sorted(grouped.items()):
        chosen = preferred_config(group)
        selected.append(chosen)
        for cfg in group:
            if cfg is not chosen:
                selected.append(replace(cfg, skip_reason=f"duplicate_disease_study_preferred:{chosen.source_disease}"))
    return selected


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
    if disease == "IBS":
        return "case" if s in IBS_CASE_STATUSES else None
    return None


def load_metadata(config: ArmConfig) -> pd.DataFrame:
    if config.metadata_path is None or not config.metadata_path.exists():
        return pd.DataFrame()
    meta = pd.read_csv(config.metadata_path, sep="\t", dtype=str)
    if "sample_alias" not in meta.columns or "subject_disease_status" not in meta.columns:
        return pd.DataFrame()
    meta = meta[["sample_alias", "subject_disease_status"]].copy()
    meta["group"] = meta["subject_disease_status"].map(lambda x: label_status(config.disease, x))
    meta = meta.dropna(subset=["group"])
    return meta.drop_duplicates("sample_alias", keep="first")


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


def base_inventory(config: ArmConfig) -> dict[str, object]:
    def release_path(path: Path | None) -> str:
        if path is None:
            return ""
        try:
            return str(path.resolve().relative_to(ROOT.resolve()))
        except ValueError:
            return str(path)

    return {
        "disease": config.disease,
        "source_disease": config.source_disease,
        "study_code": config.study_code,
        "arm_id": config.arm_id,
        "comparison": config.comparison,
        "metadata_path": release_path(config.metadata_path),
        "profile_path": release_path(config.profile_path),
        "metadata_status_counts": metadata_counts(config.metadata_path),
        "metadata_rows_all": 0,
        "labeled_samples": 0,
        "profile_samples": 0,
        "profile_species": 0,
        "matched_samples": 0,
        "case_samples": 0,
        "control_samples": 0,
        "filtered_species": 0,
        "included_main": False,
        "exclusion_reason": config.skip_reason,
        "cv_auc_mean": np.nan,
        "cv_auc_sd": np.nan,
    }


def analyze_arm(config: ArmConfig, n_estimators: int, n_jobs: int) -> tuple[pd.DataFrame, dict[str, object]]:
    inv = base_inventory(config)
    if config.skip_reason:
        return pd.DataFrame(), inv

    if config.metadata_path is None or not config.metadata_path.exists():
        inv["exclusion_reason"] = "missing_metadata_all_wide"
        return pd.DataFrame(), inv
    try:
        inv["metadata_rows_all"] = sum(1 for _ in config.metadata_path.open()) - 1
    except OSError:
        inv["metadata_rows_all"] = 0

    if config.profile_path is None or not config.profile_path.exists():
        inv["exclusion_reason"] = "missing_metaphlan4_species_profile"
        return pd.DataFrame(), inv

    meta = load_metadata(config)
    inv["labeled_samples"] = int(meta.shape[0])
    if meta.empty:
        inv["exclusion_reason"] = "missing_or_unusable_metadata_labels"
        return pd.DataFrame(), inv

    X = load_metaphlan_species(config.profile_path)
    inv["profile_samples"] = int(X.shape[0])
    inv["profile_species"] = int(X.shape[1])

    joined = meta.set_index("sample_alias").join(X, how="inner")
    inv["matched_samples"] = int(joined.shape[0])
    if joined.empty:
        inv["exclusion_reason"] = "no_profile_metadata_overlap"
        return pd.DataFrame(), inv

    y = (joined["group"] == "case").astype(int).to_numpy()
    X = joined.drop(columns=["subject_disease_status", "group"]).astype(np.float32)
    case_n = int((y == 1).sum())
    ctrl_n = int((y == 0).sum())
    inv["case_samples"] = case_n
    inv["control_samples"] = ctrl_n
    if case_n < MIN_CASE or ctrl_n < MIN_CONTROL:
        inv["exclusion_reason"] = f"insufficient_case_or_control_samples(case={case_n},control={ctrl_n})"
        return pd.DataFrame(), inv

    detection_n = (X > 0).sum(axis=0)
    mean_abund = X.mean(axis=0)
    keep = ((detection_n / X.shape[0] >= MIN_DETECTION_FRAC) | (detection_n >= MIN_DETECTION_N)) & (
        mean_abund >= MIN_MEAN_ABUND
    )
    Xf = X.loc[:, keep]
    inv["filtered_species"] = int(Xf.shape[1])
    if Xf.shape[1] == 0:
        inv["exclusion_reason"] = "no_species_after_filtering"
        return pd.DataFrame(), inv

    case_X = Xf.iloc[y == 1, :]
    ctrl_X = Xf.iloc[y == 0, :]
    species = np.array(Xf.columns)
    pvals = np.empty(len(species), dtype=float)
    case_mean = case_X.mean(axis=0).to_numpy(dtype=float)
    ctrl_mean = ctrl_X.mean(axis=0).to_numpy(dtype=float)
    case_prev = (case_X > 0).mean(axis=0).to_numpy(dtype=float)
    ctrl_prev = (ctrl_X > 0).mean(axis=0).to_numpy(dtype=float)
    log2fc = np.log2((case_mean + LOG_PSEUDOCOUNT) / (ctrl_mean + LOG_PSEUDOCOUNT))

    for i, sp in enumerate(species):
        try:
            pvals[i] = mannwhitneyu(case_X[sp].to_numpy(), ctrl_X[sp].to_numpy(), alternative="two-sided").pvalue
        except ValueError:
            pvals[i] = 1.0
    qvals = bh_fdr(pvals)

    n_splits = min(5, case_n, ctrl_n)
    aucs: list[float] = []
    importances = np.zeros(len(species), dtype=float)
    top50_counts = np.zeros(len(species), dtype=float)
    if n_splits >= 2:
        skf = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=RANDOM_SEED)
        for fold_idx, (train_idx, test_idx) in enumerate(skf.split(Xf, y)):
            clf = RandomForestClassifier(
                n_estimators=n_estimators,
                class_weight="balanced_subsample",
                max_features="sqrt",
                min_samples_leaf=1,
                random_state=RANDOM_SEED + fold_idx,
                n_jobs=n_jobs,
            )
            clf.fit(Xf.iloc[train_idx, :], y[train_idx])
            probs = clf.predict_proba(Xf.iloc[test_idx, :])[:, 1]
            if len(np.unique(y[test_idx])) == 2:
                aucs.append(float(roc_auc_score(y[test_idx], probs)))
            fold_importance = clf.feature_importances_
            importances += fold_importance / n_splits
            fold_rank = pd.Series(fold_importance, index=species).rank(ascending=False, method="min").to_numpy()
            top50_counts += (fold_rank <= min(TOP_N_MAIN, len(species))) / n_splits

    rank = pd.Series(importances, index=species).rank(ascending=False, method="min").to_numpy(dtype=int)
    ml_top_cutoff = max(TOP_N_MAIN, int(math.ceil(0.05 * len(species))))

    inv["included_main"] = True
    inv["exclusion_reason"] = ""
    inv["cv_auc_mean"] = float(np.mean(aucs)) if aucs else np.nan
    inv["cv_auc_sd"] = float(np.std(aucs, ddof=1)) if len(aucs) > 1 else np.nan

    rows = pd.DataFrame(
        {
            "disease": config.disease,
            "source_disease": config.source_disease,
            "study_code": config.study_code,
            "arm_id": config.arm_id,
            "species": species,
            "case_samples": case_n,
            "control_samples": ctrl_n,
            "n_samples": case_n + ctrl_n,
            "case_mean": case_mean,
            "control_mean": ctrl_mean,
            "case_prevalence": case_prev,
            "control_prevalence": ctrl_prev,
            "log2fc": log2fc,
            "p_value": pvals,
            "q_value": qvals,
            "rf_importance": importances,
            "rf_rank": rank,
            "rf_top50_stability": top50_counts,
            "rf_top_cutoff": ml_top_cutoff,
            "is_ml_top": rank <= ml_top_cutoff,
        }
    )
    rows["is_diff_support"] = (rows["q_value"] <= DIFF_Q_CUTOFF) & (rows["log2fc"].abs() >= LOG2FC_CUTOFF)
    rows["is_arm_key_strict"] = rows["is_diff_support"] & rows["is_ml_top"]
    return rows, inv


def aggregate_disease_keys(arm_rows: pd.DataFrame, included_inventory: pd.DataFrame) -> pd.DataFrame:
    out_rows: list[dict[str, object]] = []
    if arm_rows.empty:
        return pd.DataFrame()

    for disease, sub in arm_rows.groupby("disease"):
        n_total = int(included_inventory.loc[included_inventory["disease"] == disease, "arm_id"].nunique())
        required_support = 1 if disease == "IBS" else 2
        for species, g in sub.groupby("species"):
            weights = np.sqrt(g["n_samples"].to_numpy(dtype=float))
            p = np.clip(g["p_value"].to_numpy(dtype=float), 1e-300, 1.0)
            signed_z = np.sign(g["log2fc"].to_numpy(dtype=float)) * norm.isf(p / 2.0)
            weight_denom = math.sqrt(float(np.sum(weights**2)))
            meta_z = float(np.sum(weights * signed_z) / weight_denom) if weight_denom else 0.0
            meta_p = float(2 * norm.sf(abs(meta_z)))
            meta_lfc = float(np.average(g["log2fc"], weights=weights)) if weights.sum() else float(g["log2fc"].mean())
            signs = np.sign(g["log2fc"].to_numpy(dtype=float))
            pos = int((signs > 0).sum())
            neg = int((signs < 0).sum())
            direction_consistency = max(pos, neg) / max(pos + neg, 1)
            ml_stability = float(g["rf_top50_stability"].mean())
            n_ml_support = int((g["rf_top50_stability"] > 0).sum())
            out_rows.append(
                {
                    "disease": disease,
                    "species": species,
                    "n_arms_tested": int(g["arm_id"].nunique()),
                    "n_arms_total": n_total,
                    "n_diff_support": int(g["is_diff_support"].sum()),
                    "n_ml_support": n_ml_support,
                    "required_support": required_support,
                    "direction": "up" if meta_lfc > 0 else "down",
                    "direction_consistency": direction_consistency,
                    "meta_log2fc": meta_lfc,
                    "meta_z": meta_z,
                    "meta_p": meta_p,
                    "median_rf_rank": float(g["rf_rank"].median()),
                    "mean_rf_importance": float(g["rf_importance"].mean()),
                    "ml_stability": ml_stability,
                }
            )

    df = pd.DataFrame(out_rows)
    if df.empty:
        return df

    df["meta_q"] = np.nan
    for disease in df["disease"].unique():
        mask = df["disease"] == disease
        df.loc[mask, "meta_q"] = bh_fdr(df.loc[mask, "meta_p"].to_numpy(dtype=float))

    neglogq = -np.log10(df["meta_q"].clip(lower=1e-300))
    df["rank_score"] = neglogq * df["meta_log2fc"].abs() * df["direction_consistency"] * (1.0 + df["ml_stability"])
    df["eligible_for_top_key"] = (
        (df["direction_consistency"] >= DIRECTION_CONSISTENCY_CUTOFF)
        & (df["n_arms_tested"] >= df["required_support"])
        & ((df["disease"] != "IBS") | (df["ml_stability"] > 0))
    )
    df["ml_stability_rank"] = np.nan
    for disease in df["disease"].unique():
        mask = df["disease"] == disease
        df.loc[mask, "ml_stability_rank"] = df.loc[mask, "ml_stability"].rank(ascending=False, method="min")
    df["strict_ml_supported"] = False
    for disease, sub in df.groupby("disease"):
        cutoff = max(TOP_N_MAIN, int(math.ceil(0.10 * sub.shape[0])))
        mask = (df["disease"] == disease) & (df["ml_stability_rank"] <= cutoff) & (df["ml_stability"] > 0)
        df.loc[mask, "strict_ml_supported"] = True
    df["is_strict_key_species"] = (
        (df["meta_q"] <= DIFF_Q_CUTOFF)
        & (df["meta_log2fc"].abs() >= LOG2FC_CUTOFF)
        & (df["direction_consistency"] >= DIRECTION_CONSISTENCY_CUTOFF)
        & (df["n_diff_support"] >= df["required_support"])
        & df["strict_ml_supported"]
    )
    df["is_top50_key_species"] = False
    for disease in DISEASES:
        idx = (
            df[(df["disease"] == disease) & df["eligible_for_top_key"]]
            .sort_values("rank_score", ascending=False)
            .head(TOP_N_MAIN)
            .index
        )
        df.loc[idx, "is_top50_key_species"] = True
    return df.sort_values(["disease", "is_top50_key_species", "rank_score"], ascending=[True, False, False])


def key_sets_from_ranked(
    disease_table: pd.DataFrame,
    top_n: int = TOP_N_MAIN,
    named_only: bool = False,
    strict_only: bool = False,
) -> dict[str, set[str]]:
    key_sets: dict[str, set[str]] = {}
    for disease in DISEASES:
        sub = disease_table[disease_table["disease"] == disease].copy()
        if named_only:
            sub = sub[~sub["species"].str.contains("GGB|SGB", regex=True, na=False)]
        if strict_only:
            sub = sub[sub["is_strict_key_species"]]
        else:
            sub = sub[sub["eligible_for_top_key"]].sort_values("rank_score", ascending=False).head(top_n)
        key_sets[disease] = set(sub["species"])
    return key_sets


def background_sets(arm_rows: pd.DataFrame) -> dict[str, set[str]]:
    return {disease: set(group["species"]) for disease, group in arm_rows.groupby("disease")}


def overlap_tests(key_sets: dict[str, set[str]], backgrounds: dict[str, set[str]], permutations: int) -> pd.DataFrame:
    rng = np.random.default_rng(RANDOM_SEED)
    rows: list[dict[str, object]] = []
    for r in range(2, len(DISEASES) + 1):
        for combo in itertools.combinations(DISEASES, r):
            sets = [key_sets.get(d, set()) for d in combo]
            sizes = [len(s) for s in sets]
            shared = set.intersection(*sets) if all(sizes) else set()
            bg = set.union(*(backgrounds.get(d, set()) for d in combo))
            bg_list = np.array(sorted(bg), dtype=object)
            observed = len(shared)
            M = len(bg)
            expected_hyper = np.nan
            hyper_p = np.nan
            if r == 2 and M > 0:
                expected_hyper = sizes[0] * sizes[1] / M
                hyper_p = float(hypergeom.sf(observed - 1, M, sizes[0], sizes[1])) if observed > 0 else 1.0

            if M > 0 and all(k <= M for k in sizes):
                perm_counts = np.empty(permutations, dtype=np.int16)
                for i in range(permutations):
                    sampled = [set(rng.choice(bg_list, size=k, replace=False)) for k in sizes]
                    perm_counts[i] = len(set.intersection(*sampled))
                expected_perm = float(np.mean(perm_counts))
                perm_p = float((1 + np.sum(perm_counts >= observed)) / (permutations + 1))
            else:
                expected_perm = np.nan
                perm_p = np.nan

            union_size = len(set.union(*sets)) if sets else 0
            rows.append(
                {
                    "diseases": "|".join(combo),
                    "n_diseases": r,
                    "set_sizes": "|".join(map(str, sizes)),
                    "background_size": M,
                    "shared_count": observed,
                    "jaccard": observed / union_size if union_size else np.nan,
                    "overlap_coefficient": observed / min(sizes) if min(sizes) else np.nan,
                    "expected_hypergeom": expected_hyper,
                    "expected_permutation": expected_perm,
                    "observed_expected_ratio": observed / expected_perm if expected_perm and expected_perm > 0 else np.nan,
                    "hypergeom_p": hyper_p,
                    "permutation_p": perm_p,
                    "shared_species": ";".join(sorted(shared)),
                }
            )

    df = pd.DataFrame(rows)
    df["hypergeom_q"] = bh_fdr(df["hypergeom_p"].to_numpy(dtype=float))
    df["permutation_q"] = bh_fdr(df["permutation_p"].to_numpy(dtype=float))
    df["passes_large_overlap_rule"] = (
        (df["observed_expected_ratio"] > 2)
        & (df["permutation_q"] < 0.05)
        & (df["overlap_coefficient"] >= 0.25)
    )
    return df


def shared_species_table(disease_table: pd.DataFrame, key_sets: dict[str, set[str]]) -> pd.DataFrame:
    all_sets = [s for s in key_sets.values() if s]
    all_species = sorted(set.union(*all_sets)) if all_sets else []
    if not all_species:
        return pd.DataFrame()
    indexed = disease_table.set_index(["disease", "species"])
    rows: list[dict[str, object]] = []
    for species in all_species:
        diseases = [d for d in DISEASES if species in key_sets.get(d, set())]
        if len(diseases) < 2:
            continue
        directions: list[str] = []
        row: dict[str, object] = {
            "species": species,
            "n_diseases": len(diseases),
            "diseases": "|".join(diseases),
        }
        for disease in DISEASES:
            if disease in diseases:
                hit = indexed.loc[(disease, species)]
                directions.append(str(hit["direction"]))
                row[f"{disease}_direction"] = hit["direction"]
                row[f"{disease}_meta_log2fc"] = hit["meta_log2fc"]
                row[f"{disease}_meta_q"] = hit["meta_q"]
                row[f"{disease}_rank_score"] = hit["rank_score"]
            else:
                row[f"{disease}_direction"] = ""
                row[f"{disease}_meta_log2fc"] = np.nan
                row[f"{disease}_meta_q"] = np.nan
                row[f"{disease}_rank_score"] = np.nan
        row["direction_pattern"] = "|".join(directions)
        row["same_direction"] = len(set(directions)) == 1
        rows.append(row)
    return pd.DataFrame(rows).sort_values(["n_diseases", "same_direction", "species"], ascending=[False, False, True])


def write_plots(
    disease_table: pd.DataFrame,
    overlap: pd.DataFrame,
    shared: pd.DataFrame,
    key_sets: dict[str, set[str]],
    inventory: pd.DataFrame,
) -> None:
    FIGDIR.mkdir(parents=True, exist_ok=True)
    sns.set_theme(style="whitegrid")

    pair = overlap[overlap["n_diseases"] == 2].copy()
    mat = pd.DataFrame(np.nan, index=DISEASES, columns=DISEASES)
    for _, row in pair.iterrows():
        a, b = row["diseases"].split("|")
        mat.loc[a, b] = mat.loc[b, a] = row["overlap_coefficient"]
    np.fill_diagonal(mat.values, 1.0)
    plt.figure(figsize=(5.4, 4.3))
    sns.heatmap(mat.astype(float), annot=True, fmt=".2f", cmap="viridis", vmin=0, vmax=1)
    plt.title("Pairwise Top 50 key species overlap")
    plt.tight_layout()
    plt.savefig(FIGDIR / "pairwise_overlap_heatmap.png", dpi=220)
    plt.close()

    all_key = sorted(set.union(*(s for s in key_sets.values() if s))) if any(key_sets.values()) else []
    combo_counts = Counter("|".join(d for d in DISEASES if sp in key_sets.get(d, set())) for sp in all_key)
    if combo_counts:
        counts = pd.Series(dict(combo_counts)).sort_values(ascending=False)
        plt.figure(figsize=(max(7, len(counts) * 0.55), 4.5))
        counts.plot(kind="bar", color="#4C78A8")
        plt.ylabel("Species count")
        plt.title("Top 50 key species membership combinations")
        plt.xticks(rotation=45, ha="right")
        plt.tight_layout()
        plt.savefig(FIGDIR / "upset_style_membership_counts.png", dpi=220)
        plt.close()

    if not shared.empty:
        top_species = shared.head(40)["species"].tolist()
        rows = []
        indexed = disease_table.set_index(["disease", "species"])
        for species in top_species:
            for disease in DISEASES:
                if species not in key_sets.get(disease, set()):
                    continue
                hit = indexed.loc[(disease, species)]
                rows.append(
                    {
                        "species": species,
                        "disease": disease,
                        "meta_log2fc": hit["meta_log2fc"],
                        "neglog10q": -math.log10(max(float(hit["meta_q"]), 1e-300)),
                    }
                )
        bubble = pd.DataFrame(rows)
        if not bubble.empty:
            plt.figure(figsize=(7.2, max(4.5, len(top_species) * 0.23)))
            sns.scatterplot(
                data=bubble,
                x="disease",
                y="species",
                hue="meta_log2fc",
                size="neglog10q",
                palette="coolwarm",
                sizes=(35, 230),
                edgecolor="black",
                linewidth=0.25,
            )
            plt.title("Shared Top 50 species effect direction")
            plt.legend(bbox_to_anchor=(1.02, 1), loc="upper left", borderaxespad=0)
            plt.tight_layout()
            plt.savefig(FIGDIR / "shared_species_effect_bubble.png", dpi=220)
            plt.close()

    included = inventory[inventory["included_main"]].copy()
    if not included.empty:
        plot_df = included.sort_values(["disease", "study_code"])
        labels = plot_df["disease"] + " | " + plot_df["study_code"]
        y = np.arange(plot_df.shape[0])
        plt.figure(figsize=(9, max(5, plot_df.shape[0] * 0.28)))
        plt.barh(y, plot_df["control_samples"], color="#72B7B2", label="control")
        plt.barh(y, plot_df["case_samples"], left=plot_df["control_samples"], color="#F58518", label="case")
        plt.yticks(y, labels, fontsize=8)
        plt.xlabel("Matched samples")
        plt.title("Included disease-arm sample counts")
        plt.legend(loc="lower right")
        plt.tight_layout()
        plt.savefig(FIGDIR / "study_arm_inclusion_plot.png", dpi=220)
        plt.close()


def evaluate_full_claim(overlap: pd.DataFrame, shared: pd.DataFrame) -> tuple[bool, list[str]]:
    reasons: list[str] = []
    pair = overlap[overlap["n_diseases"] == 2]
    pass_counts = {d: 0 for d in DISEASES}
    for _, row in pair[pair["passes_large_overlap_rule"]].iterrows():
        for disease in row["diseases"].split("|"):
            pass_counts[disease] += 1
    enough_pairwise = all(pass_counts[d] >= 2 for d in DISEASES)
    if not enough_pairwise:
        reasons.append(f"pairwise_pass_counts={pass_counts}")

    multi = overlap[overlap["n_diseases"] >= 3].copy()
    has_cra_ibs_multi = bool(
        multi[
            multi["passes_large_overlap_rule"]
            & multi["diseases"].str.contains("CRA", regex=False)
            & multi["diseases"].str.contains("IBS", regex=False)
        ].shape[0]
    )
    if not has_cra_ibs_multi:
        reasons.append("no_passing_3_or_4_way_overlap_containing_CRA_and_IBS")

    shared3 = int((shared["n_diseases"] >= 3).sum()) if not shared.empty else 0
    shared4 = int((shared["n_diseases"] >= 4).sum()) if not shared.empty else 0
    enough_multi_species = shared3 >= 10 or shared4 >= 3
    if not enough_multi_species:
        reasons.append(f"multi_disease_shared_species_count_insufficient(shared3={shared3},shared4={shared4})")

    return enough_pairwise and has_cra_ibs_multi and enough_multi_species, reasons


def write_summary(
    inventory: pd.DataFrame,
    ranked: pd.DataFrame,
    top50: pd.DataFrame,
    strict: pd.DataFrame,
    overlap: pd.DataFrame,
    shared: pd.DataFrame,
    permutations: int,
) -> None:
    claim_supported, reasons = evaluate_full_claim(overlap, shared)
    lines: list[str] = []
    lines.append("# CRA/CRC/IBD/IBS Key Species Overlap V2 Report")
    lines.append("")
    lines.append("## Included disease-arms")
    included = inventory[inventory["included_main"]]
    for disease in DISEASES:
        names = included.loc[included["disease"] == disease, "study_code"].tolist()
        lines.append(f"- {disease}: {len(names)} ({', '.join(names)})")
    lines.append("")
    lines.append("## Key species counts")
    for disease in DISEASES:
        lines.append(
            f"- {disease}: Top50={int((top50['disease'] == disease).sum())}; "
            f"strict={int((strict['disease'] == disease).sum())}; "
            f"ranked={int((ranked['disease'] == disease).sum())}"
        )
    lines.append("")
    lines.append(f"Overlap tests used {permutations} permutations for the main Top50 analysis.")
    pair = overlap[overlap["n_diseases"] == 2].sort_values("permutation_q")
    for _, row in pair.iterrows():
        lines.append(
            f"- {row['diseases']}: shared={int(row['shared_count'])}, "
            f"obs/exp={row['observed_expected_ratio']:.2f}, "
            f"overlap_coeff={row['overlap_coefficient']:.2f}, "
            f"perm_q={row['permutation_q']:.4g}, "
            f"pass={bool(row['passes_large_overlap_rule'])}"
        )
    lines.append("")
    shared2 = shared.shape[0] if not shared.empty else 0
    shared3 = int((shared["n_diseases"] >= 3).sum()) if not shared.empty else 0
    shared4 = int((shared["n_diseases"] >= 4).sum()) if not shared.empty else 0
    lines.append(f"Shared Top50 species in >=2 diseases: {shared2}")
    lines.append(f"Shared Top50 species in >=3 diseases: {shared3}")
    lines.append(f"Shared Top50 species in all 4 diseases: {shared4}")
    lines.append("")
    if claim_supported:
        lines.append("Conclusion: the pre-specified criteria support the claim that CRA, CRC, IBD and IBS contain many overlapping key species.")
    else:
        lines.append("Conclusion: the pre-specified full four-disease claim is not fully supported; report the strongest passing disease-pair and multi-way results instead.")
        lines.append("Failed criteria: " + "; ".join(reasons))

    report = "\n".join(lines) + "\n"
    (OUTDIR / "analysis_summary.txt").write_text(report)
    (OUTDIR / "conclusion_report.md").write_text(report)


def run_analysis_subset(
    label: str,
    arm_rows: pd.DataFrame,
    inventory: pd.DataFrame,
    top_n: int,
    strict_only: bool,
    named_only: bool,
    permutations: int,
) -> pd.DataFrame:
    included = inventory[inventory["included_main"]].copy()
    ranked = aggregate_disease_keys(arm_rows, included)
    keys = key_sets_from_ranked(ranked, top_n=top_n, strict_only=strict_only, named_only=named_only)
    overlap = overlap_tests(keys, background_sets(arm_rows), permutations)
    overlap.insert(0, "sensitivity", label)
    return overlap


def build_sensitivity(
    arm_rows: pd.DataFrame,
    inventory: pd.DataFrame,
    ranked: pd.DataFrame,
    permutations: int,
) -> pd.DataFrame:
    rows: list[pd.DataFrame] = []
    small_perm = min(permutations, 2000)
    backgrounds = background_sets(arm_rows)

    for top_n in (25, 50, 100):
        keys = key_sets_from_ranked(ranked, top_n=top_n)
        sens = overlap_tests(keys, backgrounds, small_perm)
        sens.insert(0, "sensitivity", f"top_{top_n}")
        rows.append(sens)

    keys = key_sets_from_ranked(ranked, top_n=TOP_N_MAIN, named_only=True)
    sens = overlap_tests(keys, backgrounds, small_perm)
    sens.insert(0, "sensitivity", "top50_named_species_only")
    rows.append(sens)

    keys = key_sets_from_ranked(ranked, strict_only=True)
    sens = overlap_tests(keys, backgrounds, small_perm)
    sens.insert(0, "sensitivity", "strict_only")
    rows.append(sens)

    crc_only_rows = arm_rows[(arm_rows["disease"] != "CRC") | (arm_rows["study_code"].isin(CRC_ONLY_STUDIES))].copy()
    crc_only_inventory = inventory[
        (inventory["disease"] != "CRC") | (inventory["study_code"].isin(CRC_ONLY_STUDIES)) | (~inventory["included_main"])
    ].copy()
    rows.append(run_analysis_subset("crc_keys_crc_only_studies", crc_only_rows, crc_only_inventory, TOP_N_MAIN, False, False, small_perm))

    no_feng_rows = arm_rows[arm_rows["study_code"] != "Feng_2015_CRC_Austria"].copy()
    no_feng_inventory = inventory[(inventory["study_code"] != "Feng_2015_CRC_Austria") | (~inventory["included_main"])].copy()
    rows.append(run_analysis_subset("exclude_feng_2015_crc_austria", no_feng_rows, no_feng_inventory, TOP_N_MAIN, False, False, small_perm))

    shared_cra_crc = sorted(
        set(arm_rows.loc[arm_rows["disease"] == "CRA", "study_code"])
        & set(arm_rows.loc[arm_rows["disease"] == "CRC", "study_code"])
    )
    for study in shared_cra_crc:
        sub_rows = arm_rows[arm_rows["study_code"] != study].copy()
        sub_inventory = inventory[(inventory["study_code"] != study) | (~inventory["included_main"])].copy()
        rows.append(run_analysis_subset(f"leave_one_shared_study_out:{study}", sub_rows, sub_inventory, TOP_N_MAIN, False, False, small_perm))

    return pd.concat(rows, ignore_index=True)


def write_shared_control_table(inventory: pd.DataFrame) -> None:
    included = inventory[inventory["included_main"]].copy()
    rows: list[dict[str, object]] = []
    for study, sub in included.groupby("study_code"):
        diseases = sorted(set(sub["disease"]))
        if "CRA" in diseases and "CRC" in diseases:
            rows.append(
                {
                    "study_code": study,
                    "disease_arms": ",".join(diseases),
                    "cra_case_samples": int(sub.loc[sub["disease"] == "CRA", "case_samples"].sum()),
                    "crc_case_samples": int(sub.loc[sub["disease"] == "CRC", "case_samples"].sum()),
                    "control_samples_reused_within_study": int(sub["control_samples"].max()),
                    "note": "CRA and CRC arms are from the same Metalog study; sensitivity analyses remove shared studies and CRC keys can be restricted to CRC-only studies.",
                }
            )
    pd.DataFrame(rows).to_csv(OUTDIR / "shared_control_studies.tsv", sep="\t", index=False)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--permutations", type=int, default=10_000)
    parser.add_argument("--n-estimators", type=int, default=100)
    parser.add_argument("--n-jobs", type=int, default=4)
    args = parser.parse_args()

    OUTDIR.mkdir(parents=True, exist_ok=True)
    FIGDIR.mkdir(parents=True, exist_ok=True)
    os.environ.setdefault("MPLCONFIGDIR", str(OUTDIR / ".mplconfig"))

    arm_rows: list[pd.DataFrame] = []
    inventory_rows: list[dict[str, object]] = []
    for config in build_arm_configs():
        print(f"Analyzing {config.arm_id}", flush=True)
        rows, inventory = analyze_arm(config, args.n_estimators, args.n_jobs)
        inventory_rows.append(inventory)
        if not rows.empty:
            arm_rows.append(rows)

    inventory = pd.DataFrame(inventory_rows)
    inventory.to_csv(OUTDIR / "disease_arm_inventory.tsv", sep="\t", index=False)
    if not arm_rows:
        raise RuntimeError("No analyzable MetaPhlAn4 disease-arms were found.")

    arm_stats = pd.concat(arm_rows, ignore_index=True)
    arm_stats.to_csv(OUTDIR / "arm_level_species_stats.tsv", sep="\t", index=False)

    included = inventory[inventory["included_main"]].copy()
    ranked = aggregate_disease_keys(arm_stats, included)
    ranked.to_csv(OUTDIR / "disease_key_species_ranked.tsv", sep="\t", index=False)
    top50 = ranked[ranked["is_top50_key_species"]].copy()
    top50.to_csv(OUTDIR / "disease_key_species_top50.tsv", sep="\t", index=False)
    strict = ranked[ranked["is_strict_key_species"]].copy()
    strict.to_csv(OUTDIR / "disease_key_species_strict.tsv", sep="\t", index=False)

    key_sets = key_sets_from_ranked(ranked, top_n=TOP_N_MAIN)
    overlap = overlap_tests(key_sets, background_sets(arm_stats), args.permutations)
    overlap.to_csv(OUTDIR / "overlap_summary_top50.tsv", sep="\t", index=False)
    shared = shared_species_table(ranked, key_sets)
    shared.to_csv(OUTDIR / "shared_key_species_top50.tsv", sep="\t", index=False)

    sensitivity = build_sensitivity(arm_stats, inventory, ranked, args.permutations)
    sensitivity.to_csv(OUTDIR / "sensitivity_overlap_summary.tsv", sep="\t", index=False)
    write_shared_control_table(inventory)
    write_plots(ranked, overlap, shared, key_sets, inventory)
    write_summary(inventory, ranked, top50, strict, overlap, shared, args.permutations)

    print(f"Wrote outputs to {OUTDIR}")
    print((OUTDIR / "analysis_summary.txt").read_text())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
