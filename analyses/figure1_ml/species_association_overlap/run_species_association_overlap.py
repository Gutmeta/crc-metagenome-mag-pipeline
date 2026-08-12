#!/usr/bin/env python3
"""Species-level overlap analysis for disease-associated bacteria."""

from __future__ import annotations

import argparse
import itertools
import math
import os
from collections import Counter
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", str(Path(__file__).resolve().parent / ".mplconfig"))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy.stats import hypergeom, norm


ROOT = Path(__file__).resolve().parents[1]
OUTDIR = ROOT / "species_association_overlap"
FIGDIR = OUTDIR / "figures"
V2_DIR = ROOT / "key_species_overlap_v2"
EXTERNAL_IBS_DIR = OUTDIR / "external_ibs_associations"
DISEASES = ("CRA", "CRC", "IBD", "IBS")
PRIMARY_TOP_N = 300
RANDOM_SEED = 20260524
EXTERNAL_IBS_Q_CUTOFF = 0.10
EXTERNAL_IBS_LOG2FC_CUTOFF = 0.50


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


def load_inputs() -> tuple[pd.DataFrame, pd.DataFrame]:
    ranked = pd.read_csv(V2_DIR / "disease_key_species_ranked.tsv", sep="\t")
    arm_stats = pd.read_csv(V2_DIR / "arm_level_species_stats.tsv", sep="\t")
    ranked["eligible_for_top_key"] = ranked["eligible_for_top_key"].astype(bool)
    ranked, arm_stats = augment_ibs_with_external_associations(ranked, arm_stats)
    return ranked, arm_stats


def load_external_ibs_associations() -> pd.DataFrame:
    required = {
        "source_study",
        "species",
        "log2fc_case_vs_control",
        "p_value",
        "q_value",
        "case_samples",
        "control_samples",
    }
    frames: list[pd.DataFrame] = []
    if not EXTERNAL_IBS_DIR.exists():
        return pd.DataFrame()

    for path in sorted(EXTERNAL_IBS_DIR.glob("*.standardized.tsv")):
        if path.name.startswith("external_ibs_association_tables_combined"):
            continue
        df = pd.read_csv(path, sep="\t")
        missing = required.difference(df.columns)
        if missing:
            raise ValueError(f"{path} is missing required columns: {sorted(missing)}")
        df = df.copy()
        df["external_source_file"] = path.name
        frames.append(df)

    if not frames:
        return pd.DataFrame()

    out = pd.concat(frames, ignore_index=True, sort=False)
    out["species"] = out["species"].astype(str)
    for col in ["log2fc_case_vs_control", "p_value", "q_value", "case_samples", "control_samples"]:
        out[col] = pd.to_numeric(out[col], errors="coerce")
    out = out.dropna(subset=["species", "log2fc_case_vs_control", "p_value", "q_value"])
    out = out[out["species"].str.startswith("s__", na=False)].copy()
    out["source_pass"] = (out["q_value"] <= EXTERNAL_IBS_Q_CUTOFF) & (
        out["log2fc_case_vs_control"].abs() >= EXTERNAL_IBS_LOG2FC_CUTOFF
    )
    out["n_samples"] = out["case_samples"].fillna(0) + out["control_samples"].fillna(0)
    return out


def build_ibs_source_rows(ranked: pd.DataFrame, arm_stats: pd.DataFrame, external: pd.DataFrame) -> pd.DataFrame:
    current = ranked[ranked["disease"] == "IBS"].copy()
    ibs_arm = arm_stats[arm_stats["disease"] == "IBS"].copy()
    sample_cols = [col for col in ["case_samples", "control_samples"] if col in ibs_arm.columns]
    if sample_cols:
        sample_info = ibs_arm.groupby("species", as_index=False)[sample_cols].max()
        current = current.merge(sample_info, on="species", how="left")
    else:
        current["case_samples"] = np.nan
        current["control_samples"] = np.nan

    current_sources = pd.DataFrame(
        {
            "source_study": "Mars_2020_IBS",
            "source_type": "existing_ranked_mars_metaphlan4",
            "source_reference": "Existing key_species_overlap_v2 IBS arm",
            "comparison": "IBS vs control patient",
            "species": current["species"].astype(str),
            "log2fc_case_vs_control": pd.to_numeric(current["meta_log2fc"], errors="coerce"),
            "p_value": pd.to_numeric(current["meta_p"], errors="coerce"),
            "q_value": pd.to_numeric(current["meta_q"], errors="coerce"),
            "case_samples": pd.to_numeric(current["case_samples"], errors="coerce"),
            "control_samples": pd.to_numeric(current["control_samples"], errors="coerce"),
            "source_pass": current["eligible_for_top_key"].astype(bool),
            "external_source_file": "",
        }
    )
    current_sources["n_samples"] = current_sources["case_samples"].fillna(0) + current_sources["control_samples"].fillna(0)

    if external.empty:
        return current_sources

    external_sources = external[
        [
            "source_study",
            "source_type",
            "source_reference",
            "comparison",
            "species",
            "log2fc_case_vs_control",
            "p_value",
            "q_value",
            "case_samples",
            "control_samples",
            "source_pass",
            "external_source_file",
            "n_samples",
        ]
    ].copy()
    source_rows = pd.concat([current_sources, external_sources], ignore_index=True, sort=False)
    source_rows = source_rows.dropna(subset=["species", "log2fc_case_vs_control", "p_value"])
    source_rows["direction"] = np.where(source_rows["log2fc_case_vs_control"] > 0, "up", "down")
    source_rows["source_pass"] = source_rows["source_pass"].astype(bool)
    return source_rows


def aggregate_ibs_sources(source_rows: pd.DataFrame, ranked_columns: list[str]) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    if source_rows.empty:
        return pd.DataFrame(columns=ranked_columns)

    for species, group in source_rows.groupby("species"):
        p = np.clip(pd.to_numeric(group["p_value"], errors="coerce").fillna(1.0).to_numpy(dtype=float), 1e-300, 1.0)
        lfc = pd.to_numeric(group["log2fc_case_vs_control"], errors="coerce").fillna(0.0).to_numpy(dtype=float)
        n_samples = pd.to_numeric(group["n_samples"], errors="coerce").fillna(1.0).clip(lower=1.0).to_numpy(dtype=float)
        weights = np.sqrt(n_samples)
        signed_z = np.sign(lfc) * norm.isf(p / 2.0)
        denom = math.sqrt(float(np.sum(weights**2)))
        meta_z = float(np.sum(weights * signed_z) / denom) if denom else 0.0
        meta_p = float(2 * norm.sf(abs(meta_z)))
        meta_lfc = float(np.average(lfc, weights=weights)) if weights.sum() else float(np.mean(lfc))
        signs = np.sign(lfc)
        pos = int((signs > 0).sum())
        neg = int((signs < 0).sum())
        direction_consistency = max(pos, neg) / max(pos + neg, 1)
        source_pass = group["source_pass"].astype(bool)
        external_mask = group["source_study"].astype(str) != "Mars_2020_IBS"
        rows.append(
            {
                "disease": "IBS",
                "species": species,
                "n_arms_tested": int(group["source_study"].nunique()),
                "n_arms_total": int(source_rows["source_study"].nunique()),
                "n_diff_support": int(source_pass.sum()),
                "n_ml_support": int(source_pass[group["source_study"].astype(str).eq("Mars_2020_IBS")].sum()),
                "required_support": 1,
                "direction": "up" if meta_lfc > 0 else "down",
                "direction_consistency": direction_consistency,
                "meta_log2fc": meta_lfc,
                "meta_z": meta_z,
                "meta_p": meta_p,
                "median_rf_rank": np.nan,
                "mean_rf_importance": 0.0,
                "ml_stability": 0.0,
                "ibs_association_sources": ";".join(sorted(group["source_study"].astype(str).unique())),
                "ibs_external_source_count": int(group.loc[external_mask, "source_study"].nunique()),
                "ibs_source_pass_count": int(source_pass.sum()),
            }
        )

    out = pd.DataFrame(rows)
    out["meta_q"] = bh_fdr(out["meta_p"].to_numpy(dtype=float))
    neglogq = -np.log10(out["meta_q"].clip(lower=1e-300))
    source_bonus = 1.0 + np.log1p(out["ibs_source_pass_count"].clip(lower=0)) / math.log(4)
    out["rank_score"] = neglogq * out["meta_log2fc"].abs() * out["direction_consistency"] * source_bonus
    out["eligible_for_top_key"] = (out["direction_consistency"] >= 0.60) & (out["ibs_source_pass_count"] >= 1)
    out["ml_stability_rank"] = out["rank_score"].rank(ascending=False, method="min")
    out["strict_ml_supported"] = False
    out["is_strict_key_species"] = False
    out["is_top50_key_species"] = False
    top50_idx = out[out["eligible_for_top_key"]].sort_values("rank_score", ascending=False).head(50).index
    out.loc[top50_idx, "is_top50_key_species"] = True

    for col in ranked_columns:
        if col not in out.columns:
            out[col] = np.nan
    extra_cols = [col for col in out.columns if col not in ranked_columns]
    return out[ranked_columns + extra_cols].sort_values(
        ["is_top50_key_species", "rank_score"], ascending=[False, False]
    )


def augment_ibs_with_external_associations(
    ranked: pd.DataFrame, arm_stats: pd.DataFrame
) -> tuple[pd.DataFrame, pd.DataFrame]:
    external = load_external_ibs_associations()
    if external.empty:
        return ranked, arm_stats

    source_rows = build_ibs_source_rows(ranked, arm_stats, external)
    ibs_ranked = aggregate_ibs_sources(source_rows, list(ranked.columns))
    ranked_augmented = pd.concat([ranked[ranked["disease"] != "IBS"], ibs_ranked], ignore_index=True, sort=False)
    ranked_augmented["eligible_for_top_key"] = ranked_augmented["eligible_for_top_key"].astype(bool)

    external_background = external[["species"]].drop_duplicates().copy()
    external_background["disease"] = "IBS"
    arm_stats_augmented = pd.concat(
        [arm_stats, external_background.reindex(columns=arm_stats.columns, fill_value=np.nan)],
        ignore_index=True,
        sort=False,
    )

    source_rows.to_csv(OUTDIR / "ibs_association_source_rows_augmented.tsv", sep="\t", index=False)
    ibs_ranked.to_csv(OUTDIR / "ibs_association_augmented_ranked.tsv", sep="\t", index=False)
    summary_rows = []
    for source, group in source_rows.groupby("source_study"):
        summary_rows.append(
            {
                "source_study": source,
                "rows_total": int(group.shape[0]),
                "species_total": int(group["species"].nunique()),
                "eligible_source_pass_species": int(group.loc[group["source_pass"], "species"].nunique()),
                "up_rows": int((group["log2fc_case_vs_control"] > 0).sum()),
                "down_rows": int((group["log2fc_case_vs_control"] < 0).sum()),
                "median_abs_log2fc": float(group["log2fc_case_vs_control"].abs().median()),
            }
        )
    pd.DataFrame(summary_rows).sort_values("source_study").to_csv(
        OUTDIR / "ibs_association_external_source_summary.tsv", sep="\t", index=False
    )
    return ranked_augmented, arm_stats_augmented


def background_sets(arm_stats: pd.DataFrame) -> dict[str, set[str]]:
    return {disease: set(group["species"]) for disease, group in arm_stats.groupby("disease")}


def association_sets(
    ranked: pd.DataFrame,
    top_n: int,
    named_only: bool = False,
    require_eligible: bool = True,
) -> dict[str, set[str]]:
    sets: dict[str, set[str]] = {}
    for disease in DISEASES:
        sub = ranked[ranked["disease"] == disease].copy()
        if require_eligible:
            sub = sub[sub["eligible_for_top_key"]]
        if named_only:
            sub = sub[~sub["species"].str.contains("GGB|SGB", regex=True, na=False)]
        sub = sub.sort_values("rank_score", ascending=False).head(top_n)
        sets[disease] = set(sub["species"])
    return sets


def overlap_tests(key_sets: dict[str, set[str]], backgrounds: dict[str, set[str]], permutations: int) -> pd.DataFrame:
    rng = np.random.default_rng(RANDOM_SEED)
    rows: list[dict[str, object]] = []
    for r in range(2, len(DISEASES) + 1):
        for combo in itertools.combinations(DISEASES, r):
            sets = [key_sets.get(d, set()) for d in combo]
            sizes = [len(s) for s in sets]
            shared = set.intersection(*sets) if all(sizes) else set()
            union_keys = set.union(*sets) if sets else set()
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
                perm_q95 = float(np.quantile(perm_counts, 0.95))
            else:
                expected_perm = np.nan
                perm_p = np.nan
                perm_q95 = np.nan

            rows.append(
                {
                    "diseases": "|".join(combo),
                    "n_diseases": r,
                    "set_sizes": "|".join(map(str, sizes)),
                    "background_size": M,
                    "shared_count": observed,
                    "jaccard": observed / len(union_keys) if union_keys else np.nan,
                    "overlap_coefficient": observed / min(sizes) if min(sizes) else np.nan,
                    "expected_hypergeom": expected_hyper,
                    "expected_permutation": expected_perm,
                    "permutation_q95": perm_q95,
                    "observed_expected_ratio": observed / expected_perm if expected_perm and expected_perm > 0 else np.nan,
                    "hypergeom_p": hyper_p,
                    "permutation_p": perm_p,
                    "shared_species": ";".join(sorted(shared)),
                }
            )
    out = pd.DataFrame(rows)
    out["hypergeom_q"] = bh_fdr(out["hypergeom_p"].to_numpy(dtype=float))
    out["permutation_q"] = bh_fdr(out["permutation_p"].to_numpy(dtype=float))
    out["significantly_above_random"] = (out["permutation_q"] < 0.05) & (out["shared_count"] > out["permutation_q95"])
    return out


def shared_species_table(ranked: pd.DataFrame, key_sets: dict[str, set[str]]) -> pd.DataFrame:
    union = set.union(*(s for s in key_sets.values() if s)) if any(key_sets.values()) else set()
    idx = ranked.set_index(["disease", "species"])
    rows: list[dict[str, object]] = []
    for species in sorted(union):
        diseases = [d for d in DISEASES if species in key_sets.get(d, set())]
        if len(diseases) < 2:
            continue
        row: dict[str, object] = {
            "species": species,
            "n_diseases": len(diseases),
            "diseases": "|".join(diseases),
            "named_species": "no" if ("GGB" in species or "SGB" in species) else "yes",
        }
        directions: list[str] = []
        for disease in DISEASES:
            if disease in diseases:
                hit = idx.loc[(disease, species)]
                row[f"{disease}_direction"] = hit["direction"]
                row[f"{disease}_rank_score"] = hit["rank_score"]
                row[f"{disease}_meta_log2fc"] = hit["meta_log2fc"]
                row[f"{disease}_meta_q"] = hit["meta_q"]
                directions.append(str(hit["direction"]))
            else:
                row[f"{disease}_direction"] = ""
                row[f"{disease}_rank_score"] = np.nan
                row[f"{disease}_meta_log2fc"] = np.nan
                row[f"{disease}_meta_q"] = np.nan
        row["direction_pattern"] = "|".join(directions)
        row["same_direction"] = len(set(directions)) == 1
        rows.append(row)
    return pd.DataFrame(rows).sort_values(["n_diseases", "same_direction", "species"], ascending=[False, False, True])


def direction_summary(shared: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    for diseases, sub in shared.groupby("diseases"):
        rows.append(
            {
                "diseases": diseases,
                "shared_species": int(sub.shape[0]),
                "same_direction_species": int(sub["same_direction"].sum()),
                "mixed_direction_species": int((~sub["same_direction"]).sum()),
                "named_species": int((sub["named_species"] == "yes").sum()),
            }
        )
    return pd.DataFrame(rows).sort_values(["shared_species", "diseases"], ascending=[False, True])


def topn_sensitivity(ranked: pd.DataFrame, backgrounds: dict[str, set[str]], permutations: int) -> pd.DataFrame:
    rows: list[pd.DataFrame] = []
    for top_n in (50, 100, 150, 200, 300, 500):
        sets = association_sets(ranked, top_n=top_n)
        overlap = overlap_tests(sets, backgrounds, permutations)
        overlap.insert(0, "analysis", f"top_{top_n}")
        rows.append(overlap)
    named = association_sets(ranked, top_n=PRIMARY_TOP_N, named_only=True)
    named_overlap = overlap_tests(named, backgrounds, permutations)
    named_overlap.insert(0, "analysis", f"top_{PRIMARY_TOP_N}_named_only")
    rows.append(named_overlap)
    return pd.concat(rows, ignore_index=True)


def write_set_table(key_sets: dict[str, set[str]], ranked: pd.DataFrame) -> None:
    rows: list[pd.DataFrame] = []
    for disease in DISEASES:
        sub = ranked[(ranked["disease"] == disease) & (ranked["species"].isin(key_sets[disease]))].copy()
        sub["association_set"] = f"top_{PRIMARY_TOP_N}"
        rows.append(sub.sort_values("rank_score", ascending=False))
    pd.concat(rows, ignore_index=True).to_csv(OUTDIR / "disease_associated_species_top300.tsv", sep="\t", index=False)


def write_plots(
    primary_overlap: pd.DataFrame,
    sensitivity: pd.DataFrame,
    shared: pd.DataFrame,
    key_sets: dict[str, set[str]],
) -> None:
    FIGDIR.mkdir(parents=True, exist_ok=True)
    sns.set_theme(style="whitegrid")

    pair = primary_overlap[primary_overlap["n_diseases"] == 2]
    mat = pd.DataFrame(np.nan, index=DISEASES, columns=DISEASES)
    for _, row in pair.iterrows():
        a, b = row["diseases"].split("|")
        mat.loc[a, b] = mat.loc[b, a] = row["shared_count"]
    np.fill_diagonal(mat.values, [len(key_sets[d]) for d in DISEASES])
    plt.figure(figsize=(5.4, 4.3))
    sns.heatmap(mat.astype(float), annot=True, fmt=".0f", cmap="mako")
    plt.title("Top300 associated species overlap count")
    plt.tight_layout()
    plt.savefig(FIGDIR / "pairwise_top300_overlap_counts.png", dpi=220)
    plt.close()

    curve = sensitivity[sensitivity["n_diseases"] == 2].copy()
    curve = curve[curve["analysis"].str.match(r"top_\\d+$", na=False)].copy()
    curve["top_n"] = curve["analysis"].str.replace("top_", "", regex=False).astype(int)
    plt.figure(figsize=(8, 4.8))
    sns.lineplot(data=curve, x="top_n", y="shared_count", hue="diseases", marker="o")
    plt.xlabel("Disease-associated set size")
    plt.ylabel("Shared species count")
    plt.title("Pairwise overlap increases with association set size")
    plt.tight_layout()
    plt.savefig(FIGDIR / "topn_pairwise_overlap_curve.png", dpi=220)
    plt.close()

    combo_counts = Counter("|".join(d for d in DISEASES if sp in key_sets[d]) for sp in set.union(*key_sets.values()))
    counts = pd.Series(dict(combo_counts)).sort_values(ascending=False)
    plt.figure(figsize=(max(7, len(counts) * 0.55), 4.5))
    counts.plot(kind="bar", color="#4C78A8")
    plt.ylabel("Species count")
    plt.title("Top300 associated species membership combinations")
    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()
    plt.savefig(FIGDIR / "top300_membership_counts.png", dpi=220)
    plt.close()

    if not shared.empty:
        top = shared.sort_values(["n_diseases", "same_direction"], ascending=[False, False]).head(50)
        rows = []
        for _, row in top.iterrows():
            for disease in DISEASES:
                value = row.get(f"{disease}_meta_log2fc")
                if pd.isna(value):
                    continue
                rows.append({"species": row["species"], "disease": disease, "meta_log2fc": value})
        bubble = pd.DataFrame(rows)
        plt.figure(figsize=(7.2, max(5, top.shape[0] * 0.20)))
        sns.scatterplot(
            data=bubble,
            x="disease",
            y="species",
            hue="meta_log2fc",
            palette="coolwarm",
            s=80,
            edgecolor="black",
            linewidth=0.2,
        )
        plt.title("Shared associated species effect directions")
        plt.legend(bbox_to_anchor=(1.02, 1), loc="upper left", borderaxespad=0)
        plt.tight_layout()
        plt.savefig(FIGDIR / "shared_species_direction_bubble.png", dpi=220)
        plt.close()


def write_report(
    primary_overlap: pd.DataFrame,
    shared: pd.DataFrame,
    direction: pd.DataFrame,
    sensitivity: pd.DataFrame,
    key_sets: dict[str, set[str]],
    permutations: int,
) -> None:
    pair = primary_overlap[primary_overlap["n_diseases"] == 2].copy()
    three = primary_overlap[primary_overlap["n_diseases"] == 3].copy()
    four = primary_overlap[primary_overlap["n_diseases"] == 4].copy()
    pair_sig = int(pair["significantly_above_random"].sum())
    three_nonzero = int((three["shared_count"] > 0).sum())
    three_sig = int(three["significantly_above_random"].sum())
    shared_ge3 = int((shared["n_diseases"] >= 3).sum())
    shared_four = int((shared["n_diseases"] == 4).sum())
    named_shared_ge3 = int(((shared["n_diseases"] >= 3) & (shared["named_species"] == "yes")).sum())

    conclusion_supported = pair_sig >= 3 and three_nonzero == 4 and shared_ge3 >= 10

    lines: list[str] = []
    lines.append("# Species-Level Disease Association Overlap Report")
    lines.append("")
    lines.append("## Definition")
    lines.append(
        f"Disease-associated species are the top {PRIMARY_TOP_N} eligible disease-level ranked MetaPhlAn4 species for CRA, CRC and IBD; IBS contributes all eligible ranked species if fewer than {PRIMARY_TOP_N} are available."
    )
    lines.append(f"Permutation tests used {permutations} fixed-size draws from disease-specific detectable species backgrounds.")
    external_summary_path = OUTDIR / "ibs_association_external_source_summary.tsv"
    if external_summary_path.exists():
        external_summary = pd.read_csv(external_summary_path, sep="\t")
        external_sources = [s for s in external_summary["source_study"].astype(str).tolist() if s != "Mars_2020_IBS"]
        if external_sources:
            lines.append(
                "IBS association ranking was augmented with external differential species tables from "
                + ", ".join(external_sources)
                + "."
            )
    lines.append("")
    lines.append("## Associated Set Sizes")
    for disease in DISEASES:
        lines.append(f"- {disease}: {len(key_sets[disease])}")
    lines.append("")
    lines.append("## Primary Top300 Overlap")
    for _, row in pair.iterrows():
        lines.append(
            f"- {row['diseases']}: shared={int(row['shared_count'])}, "
            f"expected={row['expected_permutation']:.2f}, "
            f"obs/exp={row['observed_expected_ratio']:.2f}, "
            f"perm_q={row['permutation_q']:.4g}, "
            f"significant={bool(row['significantly_above_random'])}"
        )
    lines.append("")
    lines.append("## Multi-Disease Sharing")
    for _, row in three.iterrows():
        lines.append(
            f"- {row['diseases']}: shared={int(row['shared_count'])}, "
            f"obs/exp={row['observed_expected_ratio']:.2f}, perm_q={row['permutation_q']:.4g}"
        )
    if not four.empty:
        row = four.iloc[0]
        lines.append(
            f"- {row['diseases']}: shared={int(row['shared_count'])}, "
            f"obs/exp={row['observed_expected_ratio']:.2f}, perm_q={row['permutation_q']:.4g}"
        )
    lines.append("")
    lines.append(f"Species shared by >=2 diseases: {shared.shape[0]}")
    lines.append(f"Species shared by >=3 diseases: {shared_ge3}")
    lines.append(f"Named species shared by >=3 diseases: {named_shared_ge3}")
    lines.append(f"Species shared by all 4 diseases: {shared_four}")
    lines.append("")
    lines.append("## Direction Summary")
    for _, row in direction.head(12).iterrows():
        lines.append(
            f"- {row['diseases']}: shared={int(row['shared_species'])}, "
            f"same_direction={int(row['same_direction_species'])}, "
            f"mixed_direction={int(row['mixed_direction_species'])}, named={int(row['named_species'])}"
        )
    lines.append("")
    if conclusion_supported:
        lines.append(
            "Conclusion: The data support broad species-level overlap among bacteria associated with CRA, CRC, IBD and IBS, mainly through significant pairwise overlap and nonzero three-disease sharing. Four-disease shared species exist but are limited in number."
        )
    else:
        lines.append(
            "Conclusion: The data support partial species-level overlap, but not the stronger pre-specified broad-overlap criterion across all four diseases."
        )
    lines.append("")
    lines.append("Important limitation: this analysis supports overlap among disease-associated species, not overlap among strict key species.")

    report = "\n".join(lines) + "\n"
    (OUTDIR / "conclusion_report.md").write_text(report)
    (OUTDIR / "analysis_summary.txt").write_text(report)

    compact = sensitivity[
        sensitivity["analysis"].isin(["top_100", "top_200", "top_300", "top_500", f"top_{PRIMARY_TOP_N}_named_only"])
    ].copy()
    compact.to_csv(OUTDIR / "selected_sensitivity_overlap_summary.tsv", sep="\t", index=False)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--permutations", type=int, default=10_000)
    args = parser.parse_args()

    OUTDIR.mkdir(parents=True, exist_ok=True)
    FIGDIR.mkdir(parents=True, exist_ok=True)

    ranked, arm_stats = load_inputs()
    backgrounds = background_sets(arm_stats)
    key_sets = association_sets(ranked, PRIMARY_TOP_N)

    write_set_table(key_sets, ranked)
    primary_overlap = overlap_tests(key_sets, backgrounds, args.permutations)
    primary_overlap.to_csv(OUTDIR / "overlap_summary_top300.tsv", sep="\t", index=False)
    shared = shared_species_table(ranked, key_sets)
    shared.to_csv(OUTDIR / "shared_associated_species_top300.tsv", sep="\t", index=False)
    direction = direction_summary(shared)
    direction.to_csv(OUTDIR / "direction_summary_top300.tsv", sep="\t", index=False)
    sensitivity = topn_sensitivity(ranked, backgrounds, min(args.permutations, 5000))
    sensitivity.to_csv(OUTDIR / "topn_sensitivity_overlap_summary.tsv", sep="\t", index=False)

    write_plots(primary_overlap, sensitivity, shared, key_sets)
    write_report(primary_overlap, shared, direction, sensitivity, key_sets, args.permutations)

    print(f"Wrote outputs to {OUTDIR}")
    print((OUTDIR / "analysis_summary.txt").read_text())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
