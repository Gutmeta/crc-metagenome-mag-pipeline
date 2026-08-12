#!/usr/bin/env python3
"""Original-style AUC/overlap figure with CRA and CRC merged for overlap panels."""

from __future__ import annotations

import importlib.util
import itertools
import math
import os
import shutil
from collections import Counter
from pathlib import Path

os.environ.setdefault(
    "MPLCONFIGDIR",
    str(Path(__file__).resolve().parents[1] / ".mplconfig"),
)

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
import numpy as np
import pandas as pd
from scipy.stats import hypergeom


ROOT = Path(__file__).resolve().parents[3]
REQ = ROOT / "requested_topjournal_figures"
ASSOC = ROOT / "species_association_overlap"
KEY = ROOT / "key_species_overlap_v2"
OUT = ROOT / "topjournal_robustness_v1_20260530" / "crc_cra_merged_overlap"
ROBUST_TABLES = ROOT / "topjournal_robustness_v1_20260530" / "tables"
TABLES = OUT / "tables"
FIGURES = OUT / "figures"
REPORTS = OUT / "reports"
SCRIPTS = OUT / "scripts"

ORIGINAL_SCRIPT = REQ / "plot_feature_panel_auc_overlap_story.py"

OVERLAP_SETS = ("CRC/CRA", "IBD", "IBS")
SOURCE_DISEASES = ("CRA", "CRC", "IBD", "IBS")
PANEL_C_ORDER = ("CRC/CRA|IBD", "CRC/CRA|IBS", "IBD|IBS", "CRC/CRA|IBD|IBS")
RANDOM_SEED = 20260530
N_PERMUTATIONS = 10_000

COLORS = {
    "CRC/CRA": "#0B3C68",
    "IBD": "#8ABBD3",
    "IBS": "#E9BFAE",
}
PAIRWISE_CMAP = LinearSegmentedColormap.from_list(
    "merged_overlap_blue_peach",
    ["#F7FAFC", "#E9F0F4", "#BFD9E6", "#6BA7C6", "#0B3C68"],
)
TEXT = "#111820"
MUTED = "#52616E"
GRID = "#D9E1E8"
EXPECTED = "#D7A187"
SIG_BAR = "#0B3C68"
NONSIG_BAR = "#BFD9E6"


def load_original_module():
    spec = importlib.util.spec_from_file_location("feature_panel_original", ORIGINAL_SCRIPT)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot import {ORIGINAL_SCRIPT}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def ensure_dirs() -> None:
    for path in (OUT, TABLES, FIGURES, REPORTS, SCRIPTS):
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


def load_auc_inputs() -> tuple[pd.DataFrame, pd.DataFrame]:
    auc = pd.read_csv(REQ / "figure1_final_signature_auc_by_cohort.tsv", sep="\t")
    auc = auc[auc["status"].eq("ok")].copy()
    auc["disease"] = pd.Categorical(auc["disease"], categories=("CRA", "CRC", "IBD"), ordered=True)
    auc = auc.sort_values(["disease", "auc", "left_out_cohort"], ascending=[True, True, True]).reset_index(drop=True)
    pred = pd.read_csv(REQ / "figure1_final_signature_predictions.tsv", sep="\t")
    pred = pred[pred["disease"].isin(("CRA", "CRC", "IBD"))].copy()
    pred["predicted_probability"] = pd.to_numeric(pred["predicted_probability"], errors="coerce")
    pred = pred.dropna(subset=["predicted_probability"])
    return auc, pred


def load_auc_pairwise_tests() -> pd.DataFrame:
    return pd.read_csv(ROBUST_TABLES / "auc_disease_pairwise_tests.tsv", sep="\t")


def auc_q_label(q: float) -> str:
    if not np.isfinite(q):
        return "q=NA"
    if q < 0.001:
        return "q<0.001"
    if q < 0.01:
        return f"q={q:.4f}"
    return f"q={q:.2f}"


def update_panel_a_labels(fig: plt.Figure, auc_tests: pd.DataFrame) -> None:
    for ax in fig.axes:
        for text in ax.texts:
            text.set_text(text.get_text().replace("best K=", "final K="))

    q_by_comparison = {
        str(row["comparison"]): float(row["exact_permutation_q"])
        for _, row in auc_tests.iterrows()
    }
    line = (
        "Cohort-level AUC comparison: "
        f"CRA<CRC {auc_q_label(q_by_comparison.get('CRA_vs_CRC', np.nan))}; "
        f"CRA<IBD {auc_q_label(q_by_comparison.get('CRA_vs_IBD', np.nan))}; "
        f"CRC vs IBD {auc_q_label(q_by_comparison.get('CRC_vs_IBD', np.nan))}"
    )
    fig.text(
        0.126,
        0.974,
        line,
        ha="left",
        va="top",
        fontsize=8.5,
        color=MUTED,
    )


def load_external_ibs_species() -> set[str]:
    external_dir = ASSOC / "external_ibs_associations"
    if not external_dir.exists():
        return set()
    species: set[str] = set()
    for path in sorted(external_dir.glob("*.standardized.tsv")):
        if path.name.startswith("external_ibs_association_tables_combined"):
            continue
        df = pd.read_csv(path, sep="\t", usecols=["species"])
        species.update(df["species"].dropna().astype(str))
    return {sp for sp in species if sp.startswith("s__")}


def load_sets_and_backgrounds() -> tuple[dict[str, set[str]], dict[str, set[str]], pd.DataFrame]:
    assoc = pd.read_csv(ASSOC / "disease_associated_species_top300.tsv", sep="\t")
    original_sets = {
        disease: set(assoc.loc[assoc["disease"].eq(disease), "species"].astype(str))
        for disease in SOURCE_DISEASES
    }
    merged_sets = {
        "CRC/CRA": original_sets["CRA"] | original_sets["CRC"],
        "IBD": original_sets["IBD"],
        "IBS": original_sets["IBS"],
    }

    arm_stats = pd.read_csv(KEY / "arm_level_species_stats.tsv", sep="\t", usecols=["disease", "species"])
    external_species = load_external_ibs_species()
    if external_species:
        external_bg = pd.DataFrame({"disease": "IBS", "species": sorted(external_species)})
        arm_stats = pd.concat([arm_stats, external_bg], ignore_index=True)
    original_backgrounds = {
        disease: set(group["species"].dropna().astype(str))
        for disease, group in arm_stats.groupby("disease")
    }
    merged_backgrounds = {
        "CRC/CRA": original_backgrounds.get("CRA", set()) | original_backgrounds.get("CRC", set()),
        "IBD": original_backgrounds.get("IBD", set()),
        "IBS": original_backgrounds.get("IBS", set()),
    }
    return merged_sets, merged_backgrounds, assoc


def membership_counts(sets: dict[str, set[str]]) -> pd.DataFrame:
    counter: Counter[tuple[str, ...]] = Counter()
    for species in set().union(*sets.values()):
        combo = tuple(disease for disease in OVERLAP_SETS if species in sets[disease])
        counter[combo] += 1
    rows = [
        {"diseases": "|".join(combo), "n_diseases": len(combo), "species_count": count}
        for combo, count in counter.items()
    ]
    return (
        pd.DataFrame(rows)
        .sort_values(["n_diseases", "species_count", "diseases"], ascending=[False, False, True])
        .reset_index(drop=True)
    )


def overlap_tests(sets: dict[str, set[str]], backgrounds: dict[str, set[str]]) -> pd.DataFrame:
    rng = np.random.default_rng(RANDOM_SEED)
    rows = []
    for r in range(2, len(OVERLAP_SETS) + 1):
        for combo in itertools.combinations(OVERLAP_SETS, r):
            combo_sets = [sets[d] for d in combo]
            sizes = [len(s) for s in combo_sets]
            shared = set.intersection(*combo_sets) if all(sizes) else set()
            union_keys = set.union(*combo_sets) if combo_sets else set()
            bg = set.union(*(backgrounds.get(d, set()) for d in combo))
            bg_list = np.array(sorted(bg), dtype=object)
            observed = len(shared)
            expected_hyper = np.nan
            hyper_p = np.nan
            if r == 2 and len(bg) > 0:
                expected_hyper = sizes[0] * sizes[1] / len(bg)
                hyper_p = float(hypergeom.sf(observed - 1, len(bg), sizes[0], sizes[1])) if observed > 0 else 1.0
            if len(bg_list) > 0 and all(size <= len(bg_list) for size in sizes):
                perm_counts = np.empty(N_PERMUTATIONS, dtype=np.int16)
                for i in range(N_PERMUTATIONS):
                    sampled = [set(rng.choice(bg_list, size=size, replace=False)) for size in sizes]
                    perm_counts[i] = len(set.intersection(*sampled))
                expected_perm = float(np.mean(perm_counts))
                perm_p = float((1 + np.sum(perm_counts >= observed)) / (N_PERMUTATIONS + 1))
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
                    "background_size": len(bg),
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


def q_label(q: float) -> str:
    if not np.isfinite(q):
        return "q=NA"
    if q < 0.001:
        return "q<0.001"
    if q < 0.01:
        return f"q={q:.3f}"
    return f"q={q:.2f}"


def draw_pairwise_heatmap(ax: plt.Axes, overlap: pd.DataFrame, sets: dict[str, set[str]]) -> pd.DataFrame:
    pair_rows = overlap[overlap["n_diseases"].eq(2)].copy().reset_index(drop=True)
    ratio = pd.DataFrame(np.nan, index=OVERLAP_SETS, columns=OVERLAP_SETS)
    annot = pd.DataFrame("", index=OVERLAP_SETS, columns=OVERLAP_SETS)
    for _, row in pair_rows.iterrows():
        a, b = str(row["diseases"]).split("|")
        shared = int(row["shared_count"])
        obs_exp = float(row["observed_expected_ratio"])
        q = float(row["permutation_q"])
        mark = "*" if bool(row["significantly_above_random"]) else "n.s."
        for x, y in ((a, b), (b, a)):
            ratio.loc[x, y] = obs_exp
            annot.loc[x, y] = f"{shared}\n{obs_exp:.2g}x\n{q_label(q)} {mark}"
    for disease in OVERLAP_SETS:
        ratio.loc[disease, disease] = 1.0
        annot.loc[disease, disease] = f"{len(sets[disease])}"

    ax.imshow(
        ratio.to_numpy(dtype=float),
        cmap=PAIRWISE_CMAP,
        vmin=0.8,
        vmax=max(2.0, np.nanmax(ratio.to_numpy(dtype=float))),
        aspect="auto",
        interpolation="none",
    )
    ax.set_xticks(np.arange(len(OVERLAP_SETS)))
    ax.set_yticks(np.arange(len(OVERLAP_SETS)))
    ax.set_xticklabels(OVERLAP_SETS)
    ax.set_yticklabels(OVERLAP_SETS)
    ax.tick_params(axis="both", length=0)
    ax.set_xticks(np.arange(-0.5, len(OVERLAP_SETS), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(OVERLAP_SETS), 1), minor=True)
    ax.grid(which="minor", color="white", lw=1.3)
    ax.tick_params(which="minor", bottom=False, left=False)
    for i, disease_i in enumerate(OVERLAP_SETS):
        for j, disease_j in enumerate(OVERLAP_SETS):
            color = "white" if ratio.loc[disease_i, disease_j] > 1.45 else TEXT
            ax.text(
                j,
                i,
                annot.loc[disease_i, disease_j],
                ha="center",
                va="center",
                fontsize=8.2,
                color=color,
                linespacing=1.02,
            )
    ax.set_title("Pairwise enrichment after merging CRC/CRA", loc="left", fontweight="normal", color=TEXT, pad=6)
    ax.set_box_aspect(0.68)
    ax.set_anchor("W")
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.text(-0.18, 1.08, "B", transform=ax.transAxes, ha="left", va="top", fontsize=17.0, color=TEXT)
    return pair_rows


def draw_panel_c(ax: plt.Axes, overlap: pd.DataFrame) -> pd.DataFrame:
    panel_c = overlap[overlap["diseases"].isin(PANEL_C_ORDER)].copy()
    panel_c["panel_order"] = panel_c["diseases"].map({name: i for i, name in enumerate(PANEL_C_ORDER)})
    panel_c = panel_c.sort_values("panel_order").reset_index(drop=True)
    x = np.array([0.0, 1.45, 2.9, 4.6], dtype=float)
    colors = [SIG_BAR if bool(sig) else NONSIG_BAR for sig in panel_c["significantly_above_random"]]
    ax.bar(
        x,
        panel_c["shared_count"],
        color=colors,
        edgecolor="white",
        linewidth=0.85,
        width=0.62,
        label="Observed",
        zorder=3,
    )
    ax.scatter(
        x,
        panel_c["expected_permutation"],
        color=EXPECTED,
        s=46,
        edgecolor="white",
        linewidth=0.55,
        zorder=4,
        label="Expected",
    )
    ymax = max(float(panel_c["shared_count"].max()), float(panel_c["expected_permutation"].max())) if not panel_c.empty else 1.0
    for xi, (_, row) in zip(x, panel_c.iterrows()):
        mark = "*" if bool(row["significantly_above_random"]) else "n.s."
        label = f"{float(row['observed_expected_ratio']):.1f}x\n{q_label(float(row['permutation_q']))} {mark}"
        ax.text(
            xi,
            float(row["shared_count"]) + ymax * 0.035,
            label,
            ha="center",
            va="bottom",
            fontsize=7.8,
            color=TEXT,
            linespacing=1.0,
        )
    ax.set_xticks(x)
    ax.set_xticklabels(
        ["CRC/CRA\nvs IBD", "CRC/CRA\nvs IBS", "IBD\nvs IBS", "CRC/CRA\nIBD\nIBS"],
        fontsize=7.6,
        linespacing=0.95,
    )
    ax.set_ylim(0, ymax * 1.28)
    ax.set_ylabel("Shared species", labelpad=4, fontsize=10.0)
    ax.set_title("Observed vs expected overlap", loc="left", fontweight="normal", color=TEXT, pad=6)
    ax.set_axisbelow(True)
    ax.grid(axis="y", color=GRID, lw=0.75, alpha=0.9)
    ax.legend(
        frameon=False,
        loc="upper right",
        bbox_to_anchor=(1.0, 1.03),
        ncol=2,
        handlelength=1.05,
        borderpad=0.1,
        columnspacing=0.8,
        fontsize=8.8,
    )
    ax.set_box_aspect(0.68)
    ax.set_anchor("W")
    ax.text(-0.22, 1.16, "C", transform=ax.transAxes, ha="left", va="top", fontsize=17.0, color=TEXT)
    ax.spines[["top", "right"]].set_visible(False)
    return panel_c.drop(columns=["panel_order"])


def draw_overlap_panel(
    fig: plt.Figure, spec: GridSpecFromSubplotSpec, overlap: pd.DataFrame, sets: dict[str, set[str]]
) -> tuple[pd.DataFrame, pd.DataFrame]:
    panel = GridSpecFromSubplotSpec(2, 1, subplot_spec=spec, height_ratios=[1.0, 1.0], hspace=0.30)
    ax_heatmap = fig.add_subplot(panel[0])
    ax_panel_c = fig.add_subplot(panel[1])
    pairwise = draw_pairwise_heatmap(ax_heatmap, overlap, sets)
    panel_c = draw_panel_c(ax_panel_c, overlap)
    return pairwise, panel_c


def write_reports(
    auc_plot: pd.DataFrame,
    counts: pd.DataFrame,
    overlap: pd.DataFrame,
    sets: dict[str, set[str]],
    auc_tests: pd.DataFrame,
) -> None:
    med = auc_plot.groupby("disease", observed=True)["auc"].median().to_dict()
    final_k = auc_plot.groupby("disease", observed=True)["final_k"].first().to_dict()
    pair = overlap[overlap["n_diseases"].eq(2)].copy()
    three = overlap[overlap["n_diseases"].eq(3)].iloc[0]
    lines = [
        "# CRC/CRA-Merged AUC and Overlap Figure",
        "",
        "Panel A keeps the locked CRA, CRC and IBD prediction panels unchanged, but labels the locked sizes as final K and displays cohort-level AUC comparison q values.",
        f"Final panel sizes are CRA k={int(final_k['CRA'])}, CRC k={int(final_k['CRC'])}, IBD k={int(final_k['IBD'])}; median AUCs are CRA={med['CRA']:.3f}, CRC={med['CRC']:.3f}, IBD={med['IBD']:.3f}.",
        "AUC comparison q values are BH-adjusted exact permutation q values from the topjournal robustness run.",
    ]
    for _, row in auc_tests.iterrows():
        lines.append(
            f"- {row['comparison']}: mean AUC difference={float(row['mean_diff_b_minus_a']):.3f}, exact_permutation_q={float(row['exact_permutation_q']):.4g}."
        )
    lines += [
        "",
        "Panels B-C merge CRA and CRC disease-associated Top300 species into one CRC/CRA state using the union of the original CRA and CRC associated sets.",
        f"Set sizes: CRC/CRA={len(sets['CRC/CRA'])}, IBD={len(sets['IBD'])}, IBS={len(sets['IBS'])}.",
        "",
        "## Pairwise overlap",
    ]
    for _, row in pair.iterrows():
        lines.append(
            f"- {row['diseases']}: shared={int(row['shared_count'])}, obs/exp={float(row['observed_expected_ratio']):.2f}, permutation_q={float(row['permutation_q']):.4g}, significant={bool(row['significantly_above_random'])}."
        )
    lines += [
        "",
        "## Three-way overlap",
        f"- {three['diseases']}: shared={int(three['shared_count'])}, obs/exp={float(three['observed_expected_ratio']):.2f}, permutation_q={float(three['permutation_q']):.4g}, significant={bool(three['significantly_above_random'])}.",
        "",
        "Panel C now displays the three pairwise comparisons plus the three-disease intersection together as observed-vs-expected overlap counts.",
        "",
        "Important: the overlap panels use broader disease-associated ranked species sets, not the locked prediction panels from Panel A.",
    ]
    (REPORTS / "crc_cra_merged_overlap_legend.md").write_text("\n".join(lines) + "\n")


def main() -> int:
    ensure_dirs()
    source_script = Path(__file__).resolve()
    snapshot_path = SCRIPTS / source_script.name
    if source_script != snapshot_path:
        shutil.copy2(source_script, snapshot_path)
    original = load_original_module()
    original.setup_style()

    auc, predictions = load_auc_inputs()
    auc_tests = load_auc_pairwise_tests()
    sets, backgrounds, assoc = load_sets_and_backgrounds()
    counts = membership_counts(sets)
    overlap = overlap_tests(sets, backgrounds)

    fig = plt.figure(figsize=(18.0, 5.0))
    outer = GridSpec(1, 2, figure=fig, width_ratios=[2.43, 0.94], wspace=0.085)
    auc_plot = original.draw_auc_forest(fig, outer[0], auc, predictions)
    update_panel_a_labels(fig, auc_tests)
    pairwise, panel_c = draw_overlap_panel(fig, outer[1], overlap, sets)
    fig.subplots_adjust(left=0.062, right=0.986, top=0.885, bottom=0.17)

    stem = "feature_panel_auc_overlap_story_crc_cra_merged"
    fig.savefig(FIGURES / f"{stem}.pdf")
    fig.savefig(FIGURES / f"{stem}.svg")
    fig.savefig(FIGURES / f"{stem}.png", dpi=300)
    fig.savefig(FIGURES / f"{stem}.tiff", dpi=600, pil_kwargs={"compression": "tiff_lzw"})
    plt.close(fig)

    auc_plot.drop(columns=["plot_y"], errors="ignore").to_csv(TABLES / f"{stem}_auc_source.tsv", sep="\t", index=False)
    auc_tests.to_csv(TABLES / f"{stem}_auc_pairwise_tests.tsv", sep="\t", index=False)
    counts.to_csv(TABLES / f"{stem}_membership_counts.tsv", sep="\t", index=False)
    overlap.to_csv(TABLES / f"{stem}_overlap_summary.tsv", sep="\t", index=False)
    pairwise.to_csv(TABLES / f"{stem}_pairwise_overlap.tsv", sep="\t", index=False)
    panel_c.to_csv(TABLES / f"{stem}_panel_c_overlap.tsv", sep="\t", index=False)
    overlap[overlap["n_diseases"].eq(3)].to_csv(TABLES / f"{stem}_threeway_overlap.tsv", sep="\t", index=False)
    pd.DataFrame(
        [{"merged_set": name, "species_count": len(species), "species": ";".join(sorted(species))} for name, species in sets.items()]
    ).to_csv(TABLES / f"{stem}_set_membership.tsv", sep="\t", index=False)
    write_reports(auc_plot, counts, overlap, sets, auc_tests)
    print(f"Wrote CRC/CRA-merged figure package to {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
