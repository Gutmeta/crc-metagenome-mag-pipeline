#!/usr/bin/env python3
"""Create publication-style figures for species association overlap results."""

from __future__ import annotations

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
from matplotlib import gridspec
from matplotlib.patches import FancyArrowPatch, Rectangle


ROOT = Path(__file__).resolve().parent
FIGDIR = ROOT / "publication_figures"
DISEASES = ("CRA", "CRC", "IBD", "IBS")
DISEASE_COLORS = {
    "CRA": "#4C78A8",
    "CRC": "#F58518",
    "IBD": "#54A24B",
    "IBS": "#B279A2",
}
TEXT = "#1F2933"
MUTED = "#6B7280"
GRID = "#E5E7EB"


def clean_species_name(name: str) -> str:
    return str(name).replace("s__", "").replace("_", " ")


def q_label(q: float) -> str:
    if not np.isfinite(q):
        return "n.s."
    if q < 0.001:
        return "q<0.001"
    if q < 0.01:
        return f"q={q:.3f}"
    return f"q={q:.2f}"


def setup_style() -> None:
    sns.set_theme(style="white")
    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 9,
            "axes.titlesize": 10,
            "axes.labelsize": 9,
            "xtick.labelsize": 8,
            "ytick.labelsize": 8,
            "legend.fontsize": 8,
            "axes.edgecolor": "#CBD5E1",
            "axes.linewidth": 0.8,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )


def panel_label(ax, label: str) -> None:
    ax.text(
        -0.12,
        1.08,
        label,
        transform=ax.transAxes,
        fontsize=13,
        fontweight="bold",
        va="top",
        ha="left",
        color=TEXT,
    )


def draw_framework(ax, set_sizes: dict[str, int]) -> None:
    ax.axis("off")
    panel_label(ax, "A")
    ax.text(
        0.0,
        1.02,
        "Disease-associated species were tested against a matched species background",
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=10,
        fontweight="bold",
        color=TEXT,
    )
    y = 0.70
    x_positions = [0.02, 0.25, 0.48, 0.71]
    for x, disease in zip(x_positions, DISEASES):
        ax.add_patch(
            Rectangle(
                (x, y),
                0.18,
                0.18,
                transform=ax.transAxes,
                facecolor=DISEASE_COLORS[disease],
                edgecolor="none",
                alpha=0.92,
            )
        )
        ax.text(x + 0.09, y + 0.115, disease, transform=ax.transAxes, color="white", ha="center", va="center", fontweight="bold")
        ax.text(
            x + 0.09,
            y + 0.045,
            f"Top {set_sizes[disease]}",
            transform=ax.transAxes,
            color="white",
            ha="center",
            va="center",
            fontsize=8,
        )

    ax.add_patch(FancyArrowPatch((0.17, 0.54), (0.83, 0.54), transform=ax.transAxes, arrowstyle="-|>", mutation_scale=12, lw=1.2, color=MUTED))
    ax.text(0.50, 0.59, "fixed-size permutation", transform=ax.transAxes, ha="center", color=MUTED, fontsize=8)
    ax.text(0.50, 0.49, "detectable species union as null space", transform=ax.transAxes, ha="center", color=MUTED, fontsize=8)

    summary = [
        ("196", "shared by >=2 diseases"),
        ("34", "shared by >=3 diseases"),
        ("2", "shared by all four diseases"),
    ]
    for i, (num, label) in enumerate(summary):
        x = 0.08 + i * 0.30
        ax.text(x, 0.24, num, transform=ax.transAxes, ha="center", va="center", fontsize=22, fontweight="bold", color=TEXT)
        ax.text(x, 0.13, label, transform=ax.transAxes, ha="center", va="center", fontsize=8, color=MUTED)


def draw_pairwise_enrichment(ax, overlap: pd.DataFrame) -> None:
    panel_label(ax, "B")
    pair = overlap[overlap["n_diseases"] == 2].copy()
    pair["label"] = pair["diseases"].str.replace("|", " vs ", regex=False)
    pair = pair.sort_values("observed_expected_ratio", ascending=True)
    y = np.arange(pair.shape[0])
    colors = ["#2563EB" if sig else "#9CA3AF" for sig in pair["significantly_above_random"]]
    ax.axvline(1, color="#9CA3AF", lw=1, ls="--", zorder=0)
    ax.scatter(
        pair["observed_expected_ratio"],
        y,
        s=pair["shared_count"] * 8,
        c=colors,
        edgecolor="white",
        linewidth=0.7,
        zorder=3,
    )
    for yi, (_, row) in zip(y, pair.iterrows()):
        ax.text(
            row["observed_expected_ratio"] + 0.035,
            yi,
            f"n={int(row['shared_count'])}, {q_label(row['permutation_q'])}",
            va="center",
            ha="left",
            fontsize=7.5,
            color=TEXT if row["significantly_above_random"] else MUTED,
        )
    ax.set_yticks(y)
    ax.set_yticklabels(pair["label"])
    ax.set_xlabel("Observed / expected overlap")
    ax.set_title("Five of six pairwise overlaps exceed random expectation", loc="left", fontweight="bold", color=TEXT)
    ax.set_xlim(0.75, max(1.9, pair["observed_expected_ratio"].max() + 0.45))
    ax.grid(axis="x", color=GRID, lw=0.8)
    ax.spines[["top", "right"]].set_visible(False)


def draw_multiway(ax, overlap: pd.DataFrame) -> None:
    panel_label(ax, "C")
    multi = overlap[overlap["n_diseases"] >= 3].copy()
    multi["label"] = multi["diseases"].str.replace("|", "\n", regex=False)
    multi = multi.sort_values(["n_diseases", "shared_count"], ascending=[True, False])
    x = np.arange(multi.shape[0])
    colors = ["#111827" if n == 4 else "#64748B" for n in multi["n_diseases"]]
    ax.bar(x, multi["shared_count"], color=colors, width=0.72)
    ax.scatter(x, multi["expected_permutation"], color="#F59E0B", s=24, zorder=4, label="Expected")
    for xi, (_, row) in zip(x, multi.iterrows()):
        ax.text(xi, row["shared_count"] + 0.7, q_label(row["permutation_q"]), ha="center", va="bottom", fontsize=7, color=TEXT)
    ax.set_xticks(x)
    ax.set_xticklabels(multi["label"])
    ax.set_ylabel("Shared species")
    ax.set_title("All multi-disease intersections are non-random", loc="left", fontweight="bold", color=TEXT)
    ax.grid(axis="y", color=GRID, lw=0.8)
    ax.legend(frameon=False, loc="upper right")
    ax.spines[["top", "right"]].set_visible(False)


def draw_membership(ax, shared: pd.DataFrame) -> None:
    panel_label(ax, "D")
    counts = Counter(shared["diseases"])
    ordered = pd.Series(dict(counts)).sort_values(ascending=True)
    colors = ["#111827" if label.count("|") == 3 else "#64748B" if label.count("|") == 2 else "#94A3B8" for label in ordered.index]
    ax.barh(np.arange(len(ordered)), ordered.values, color=colors)
    ax.set_yticks(np.arange(len(ordered)))
    ax.set_yticklabels(ordered.index.str.replace("|", " + ", regex=False))
    for i, value in enumerate(ordered.values):
        ax.text(value + 1, i, str(int(value)), va="center", ha="left", fontsize=8, color=TEXT)
    ax.set_xlabel("Species count")
    ax.set_title("Shared species are distributed across disease combinations", loc="left", fontweight="bold", color=TEXT)
    ax.grid(axis="x", color=GRID, lw=0.8)
    ax.set_xlim(0, max(ordered.values) + 8)
    ax.spines[["top", "right"]].set_visible(False)


def draw_direction(ax, direction: pd.DataFrame) -> None:
    panel_label(ax, "E")
    plot = direction.copy()
    plot = plot[plot["shared_species"] >= 2].sort_values("shared_species", ascending=True)
    y = np.arange(plot.shape[0])
    ax.barh(y, plot["same_direction_species"], color="#10B981", label="same direction")
    ax.barh(y, plot["mixed_direction_species"], left=plot["same_direction_species"], color="#D97706", label="mixed direction")
    ax.set_yticks(y)
    ax.set_yticklabels(plot["diseases"].str.replace("|", " + ", regex=False))
    ax.set_xlabel("Shared species")
    ax.set_title("Directionality distinguishes convergent and disease-specific signals", loc="left", fontweight="bold", color=TEXT)
    ax.legend(frameon=False, loc="lower right")
    ax.grid(axis="x", color=GRID, lw=0.8)
    ax.spines[["top", "right"]].set_visible(False)


def draw_species_heatmap(ax, shared: pd.DataFrame) -> None:
    panel_label(ax, "F")
    selected = shared[(shared["n_diseases"] >= 3) & (shared["named_species"] == "yes")].copy()
    selected["abs_score"] = selected[[f"{d}_rank_score" for d in DISEASES]].fillna(0).sum(axis=1)
    selected = selected.sort_values(["n_diseases", "same_direction", "abs_score"], ascending=[False, False, False]).head(18)
    matrix = selected.set_index("species")[[f"{d}_meta_log2fc" for d in DISEASES]]
    matrix.columns = DISEASES
    matrix.index = [clean_species_name(i) for i in matrix.index]
    sns.heatmap(
        matrix,
        ax=ax,
        cmap="vlag",
        center=0,
        vmin=-3,
        vmax=3,
        linewidths=0.4,
        linecolor="white",
        cbar_kws={"label": "meta log2FC", "shrink": 0.65},
    )
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_title("Representative multi-disease species retain effect-size structure", loc="left", fontweight="bold", color=TEXT)
    ax.tick_params(axis="y", labelsize=7)
    ax.tick_params(axis="x", rotation=0)


def make_main_figure() -> None:
    overlap = pd.read_csv(ROOT / "overlap_summary_top300.tsv", sep="\t")
    shared = pd.read_csv(ROOT / "shared_associated_species_top300.tsv", sep="\t")
    direction = pd.read_csv(ROOT / "direction_summary_top300.tsv", sep="\t")
    assoc = pd.read_csv(ROOT / "disease_associated_species_top300.tsv", sep="\t")
    set_sizes = assoc.groupby("disease")["species"].nunique().to_dict()

    setup_style()
    fig = plt.figure(figsize=(13.2, 10.4), constrained_layout=False)
    gs = gridspec.GridSpec(3, 2, figure=fig, height_ratios=[1.0, 1.25, 1.6], hspace=0.58, wspace=0.38)
    draw_framework(fig.add_subplot(gs[0, 0]), set_sizes)
    draw_pairwise_enrichment(fig.add_subplot(gs[0, 1]), overlap)
    draw_multiway(fig.add_subplot(gs[1, 0]), overlap)
    draw_membership(fig.add_subplot(gs[1, 1]), shared)
    draw_direction(fig.add_subplot(gs[2, 0]), direction)
    draw_species_heatmap(fig.add_subplot(gs[2, 1]), shared)
    fig.suptitle(
        "Disease-associated gut species converge across premalignant, malignant and inflammatory bowel phenotypes",
        x=0.02,
        y=0.995,
        ha="left",
        fontsize=14,
        fontweight="bold",
        color=TEXT,
    )
    fig.text(
        0.02,
        0.965,
        "Top-ranked MetaPhlAn4 species were evaluated against fixed-size random species backgrounds; the analysis tests association-level overlap, not strict key-species identity.",
        ha="left",
        va="top",
        fontsize=9,
        color=MUTED,
    )
    FIGDIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIGDIR / "main_species_association_overlap.pdf", bbox_inches="tight")
    fig.savefig(FIGDIR / "main_species_association_overlap.png", dpi=320, bbox_inches="tight")
    plt.close(fig)


def make_sensitivity_figure() -> None:
    sensitivity = pd.read_csv(ROOT / "topn_sensitivity_overlap_summary.tsv", sep="\t")
    setup_style()
    fig, axes = plt.subplots(1, 2, figsize=(12.5, 4.4), gridspec_kw={"width_ratios": [1.15, 1.0]})

    curve = sensitivity[(sensitivity["n_diseases"] == 2) & sensitivity["analysis"].str.match(r"top_\\d+$", na=False)].copy()
    curve["top_n"] = curve["analysis"].str.replace("top_", "", regex=False).astype(int)
    sns.lineplot(data=curve, x="top_n", y="shared_count", hue="diseases", marker="o", ax=axes[0])
    axes[0].set_title("Overlap grows smoothly with association-set breadth", loc="left", fontweight="bold", color=TEXT)
    axes[0].set_xlabel("Top ranked species per disease")
    axes[0].set_ylabel("Pairwise shared species")
    axes[0].grid(color=GRID, lw=0.8)
    axes[0].spines[["top", "right"]].set_visible(False)
    panel_label(axes[0], "A")

    named = sensitivity[sensitivity["analysis"] == "top_300_named_only"].copy()
    named = named[named["n_diseases"] >= 3].copy()
    named["label"] = named["diseases"].str.replace("|", "\n", regex=False)
    named = named.sort_values(["n_diseases", "shared_count"], ascending=[True, False])
    x = np.arange(named.shape[0])
    axes[1].bar(x, named["shared_count"], color="#334155", width=0.7)
    axes[1].scatter(x, named["expected_permutation"], color="#F59E0B", s=26, label="Expected")
    for xi, (_, row) in zip(x, named.iterrows()):
        axes[1].text(xi, row["shared_count"] + 0.4, q_label(row["permutation_q"]), ha="center", va="bottom", fontsize=7)
    axes[1].set_xticks(x)
    axes[1].set_xticklabels(named["label"])
    axes[1].set_ylabel("Named shared species")
    axes[1].set_title("Multi-disease sharing persists among named species", loc="left", fontweight="bold", color=TEXT)
    axes[1].grid(axis="y", color=GRID, lw=0.8)
    axes[1].legend(frameon=False)
    axes[1].spines[["top", "right"]].set_visible(False)
    panel_label(axes[1], "B")

    fig.suptitle("Extended Data: robustness of species-level overlap", x=0.02, y=1.03, ha="left", fontsize=13, fontweight="bold", color=TEXT)
    FIGDIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(FIGDIR / "extended_species_overlap_sensitivity.pdf", bbox_inches="tight")
    fig.savefig(FIGDIR / "extended_species_overlap_sensitivity.png", dpi=320, bbox_inches="tight")
    plt.close(fig)


def write_legend() -> None:
    legend = """# Publication Figure Legend

## Main Figure

**Disease-associated gut species converge across premalignant, malignant and inflammatory bowel phenotypes.**

**A,** Analysis schematic. Disease-associated species were defined from disease-level ranked MetaPhlAn4 profiles. CRA, CRC and IBD contributed the top 300 eligible species; IBS contributed all 81 eligible ranked species. Overlap was evaluated against fixed-size random draws from disease-specific detectable species backgrounds.

**B,** Pairwise overlap enrichment. Points show observed-to-expected overlap ratios; point area is proportional to the number of shared species. Five of six pairwise comparisons exceeded the random expectation after permutation testing.

**C,** Multi-disease sharing. Bars show observed shared species counts and orange points show expected counts under the permutation null. All three-disease intersections and the four-disease intersection were above the corresponding random expectation.

**D,** Membership distribution of shared associated species. The analysis identified 196 species shared by at least two diseases, 34 shared by at least three diseases and two shared by all four diseases.

**E,** Directional architecture of shared species. Same-direction and mixed-direction components are shown separately, indicating that shared disease association can reflect either convergent shifts or phenotype-specific directionality.

**F,** Effect-size heatmap for representative named species shared by three or more diseases. Values are disease-level meta log2 fold-changes.

## Key wording for manuscript

The strongest defensible sentence is:

> Disease-associated species-level gut bacteria showed extensive overlap across CRA, CRC, IBD and IBS, with 196 species shared by at least two diseases, 34 species shared by three or more diseases and a four-disease core of two species, exceeding fixed-size random species backgrounds.

Avoid replacing “disease-associated species” with “strict key species”; the strict key-species analysis did not support the four-disease claim.
"""
    (FIGDIR / "figure_legend.md").write_text(legend)


def main() -> int:
    make_main_figure()
    make_sensitivity_figure()
    write_legend()
    print(f"Wrote publication figures to {FIGDIR}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
