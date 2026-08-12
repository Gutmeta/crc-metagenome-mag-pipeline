#!/usr/bin/env python3
"""Create the requested publication-style figures and source tables."""

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
from matplotlib.patches import Rectangle


ROOT = Path(__file__).resolve().parents[1]
OUTDIR = ROOT / "requested_topjournal_figures"
ASSOC_DIR = ROOT / "species_association_overlap"

DISEASES = ("CRA", "CRC", "IBD", "IBS")
FIG1_DISEASES = ("CRA", "CRC", "IBD")
DISEASE_COLORS = {
    "CRA": "#3B82F6",
    "CRC": "#F97316",
    "IBD": "#16A34A",
    "IBS": "#A855F7",
}
FIG1_COLORS = {
    "CRA": "#A8C9C7",
    "CRC": "#6F7DBA",
    "IBD": "#34256F",
}
TEXT = "#111827"
MUTED = "#6B7280"
GRID = "#E5E7EB"


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


def panel_label(ax: plt.Axes, label: str) -> None:
    ax.text(
        -0.08,
        1.05,
        label,
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=13,
        fontweight="bold",
        color=TEXT,
    )


def clean_cohort(name: str) -> str:
    return str(name).replace("_CRC_", " CRC ").replace("_", " ")


def clean_species(name: str) -> str:
    return str(name).replace("s__", "").replace("_", " ")


def q_label(q: float) -> str:
    if not np.isfinite(q):
        return "n.s."
    if q < 0.001:
        return "q<0.001"
    if q < 0.01:
        return f"q={q:.3f}"
    return f"q={q:.2f}"


def p_label(p: float) -> str:
    if not np.isfinite(p):
        return "p=n.a."
    if p < 0.001:
        return "p<0.001"
    if p < 0.01:
        return f"p={p:.3f}"
    return f"p={p:.2f}"


def draw_figure1() -> None:
    source = pd.read_csv(OUTDIR / "figure1_source_data.tsv", sep="\t")
    disease_summary = pd.read_csv(OUTDIR / "figure1_disease_summary.tsv", sep="\t")
    final_auc = pd.read_csv(OUTDIR / "figure1_final_signature_auc_by_cohort.tsv", sep="\t")
    final_features = pd.read_csv(OUTDIR / "figure1_final_signature_features.tsv", sep="\t")
    final_summary = pd.read_csv(OUTDIR / "figure1_final_signature_disease_summary.tsv", sep="\t")
    k_selection = pd.read_csv(OUTDIR / "figure1_final_signature_k_selection.tsv", sep="\t")
    source = source[source["status"] == "ok"].copy()
    final_auc = final_auc[final_auc["status"] == "ok"].copy()
    if source.empty:
        raise RuntimeError("figure1_source_data.tsv has no successful AUC rows")
    if final_auc.empty:
        raise RuntimeError("figure1_final_signature_auc_by_cohort.tsv has no successful AUC rows")

    source["disease"] = pd.Categorical(source["disease"], categories=FIG1_DISEASES, ordered=True)
    source = source.sort_values(["disease", "auc", "left_out_cohort"], ascending=[True, True, True]).reset_index(drop=True)
    source["row_label"] = source.apply(lambda r: f"{r['disease']} | {clean_cohort(r['left_out_cohort'])}", axis=1)
    final_auc["disease"] = pd.Categorical(final_auc["disease"], categories=FIG1_DISEASES, ordered=True)
    final_auc = final_auc.sort_values(["disease", "auc", "left_out_cohort"], ascending=[True, True, True]).reset_index(drop=True)
    final_auc["row_label"] = final_auc.apply(lambda r: f"{r['disease']} | {clean_cohort(r['left_out_cohort'])}", axis=1)

    disease_summary["disease"] = pd.Categorical(disease_summary["disease"], categories=FIG1_DISEASES, ordered=True)
    disease_summary = disease_summary.sort_values("disease").reset_index(drop=True)
    final_summary["disease"] = pd.Categorical(final_summary["disease"], categories=FIG1_DISEASES, ordered=True)
    final_summary = final_summary.sort_values("disease").reset_index(drop=True)

    fig1_cmap = sns.blend_palette(["#E9F2F1", "#7A82BA", "#231653"], as_cmap=True)
    handles = [
        plt.Line2D([0], [0], marker="o", color="none", markerfacecolor=FIG1_COLORS[d], markeredgecolor="white", label=d, markersize=8)
        for d in FIG1_DISEASES
    ]

    fig = plt.figure(figsize=(15.2, 16.0))
    gs = gridspec.GridSpec(
        3,
        2,
        width_ratios=[3.25, 1.45],
        height_ratios=[5.0, 2.75, 5.15],
        hspace=0.36,
        wspace=0.18,
        figure=fig,
    )
    ax_auc = fig.add_subplot(gs[0, 0])
    ax_k_heat = fig.add_subplot(gs[0, 1])
    ax_summary = fig.add_subplot(gs[1, 0])
    ax_stability = fig.add_subplot(gs[1, 1])
    ax_final_auc = fig.add_subplot(gs[2, :])

    y = np.arange(source.shape[0])
    colors = [FIG1_COLORS[str(d)] for d in source["disease"]]
    for yi, (_, row) in zip(y, source.iterrows()):
        color = FIG1_COLORS[str(row["disease"])]
        ax_auc.plot([float(row["auc_ci_low"]), float(row["auc_ci_high"])], [yi, yi], color=color, lw=2.1, alpha=0.72, solid_capstyle="round", zorder=2)
        ax_auc.plot([float(row["auc_ci_low"]), float(row["auc_ci_low"])], [yi - 0.12, yi + 0.12], color=color, lw=1.3, alpha=0.72, zorder=2)
        ax_auc.plot([float(row["auc_ci_high"]), float(row["auc_ci_high"])], [yi - 0.12, yi + 0.12], color=color, lw=1.3, alpha=0.72, zorder=2)
    ax_auc.scatter(source["auc"], y, s=58, c=colors, edgecolor="white", linewidth=0.8, zorder=3)
    ax_auc.axvline(0.5, color="#9CA3AF", lw=1.0, ls="--", zorder=0)
    for yi, (_, row) in zip(y, source.iterrows()):
        sig = "*" if bool(row.get("significant_auc_gt_0_5", False)) else ""
        ax_auc.text(
            min(float(row["auc_ci_high"]) + 0.012, 1.055),
            yi,
            f"{float(row['auc']):.2f}{sig} | k={int(row['best_k'])}",
            va="center",
            ha="left",
            fontsize=7.3,
            color=TEXT,
        )
    ax_auc.set_yticks(y)
    ax_auc.set_yticklabels(source["row_label"])
    ax_auc.invert_yaxis()
    ax_auc.set_xlim(0.45, 1.10)
    ax_auc.set_xlabel("Held-out cohort AUC with bootstrap 95% CI")
    ax_auc.set_title("Nested LODO external validation with training-selected k", loc="left", fontweight="bold", color=TEXT)
    ax_auc.grid(axis="x", color=GRID, lw=0.8)
    ax_auc.spines[["top", "right"]].set_visible(False)
    panel_label(ax_auc, "A")

    k_pivot = (
        k_selection.pivot_table(index="disease", columns="candidate_k", values="mean_inner_auc", aggfunc="mean")
        .reindex(index=FIG1_DISEASES, columns=[5, 10, 15, 20, 30, 50, 75, 100, 150, 200])
    )
    sns.heatmap(
        k_pivot,
        ax=ax_k_heat,
        cmap=fig1_cmap,
        vmin=0.48,
        vmax=max(0.82, float(np.nanmax(k_pivot.to_numpy()))),
        linewidths=0.6,
        linecolor="white",
        cbar_kws={"label": "Mean inner AUC", "shrink": 0.68},
    )
    selected_k = k_selection[k_selection["selected_final_k"]].set_index("disease")["candidate_k"].to_dict()
    for row_idx, disease in enumerate(FIG1_DISEASES):
        if disease in selected_k and selected_k[disease] in k_pivot.columns:
            col_idx = list(k_pivot.columns).index(selected_k[disease])
            ax_k_heat.add_patch(Rectangle((col_idx, row_idx), 1, 1, fill=False, edgecolor="#111827", lw=1.8))
    ax_k_heat.set_xlabel("Candidate k")
    ax_k_heat.set_ylabel("")
    ax_k_heat.set_title("Best-k tuning landscape", loc="left", fontweight="bold", color=TEXT)
    panel_label(ax_k_heat, "B")

    dist_nested = source[["disease", "left_out_cohort", "auc"]].copy()
    dist_nested["evaluation"] = "Nested LODO"
    dist_locked = final_auc[["disease", "left_out_cohort", "auc"]].copy()
    dist_locked["evaluation"] = "Locked panel"
    dist = pd.concat([dist_nested, dist_locked], ignore_index=True)
    dist["disease"] = pd.Categorical(dist["disease"], categories=FIG1_DISEASES, ordered=True)
    sns.violinplot(
        data=dist,
        x="disease",
        y="auc",
        hue="evaluation",
        split=True,
        inner=None,
        cut=0,
        linewidth=0.7,
        palette={"Nested LODO": "#C9DCDA", "Locked panel": "#34256F"},
        ax=ax_summary,
    )
    sns.stripplot(
        data=dist,
        x="disease",
        y="auc",
        hue="evaluation",
        dodge=True,
        size=3.4,
        alpha=0.85,
        linewidth=0.4,
        edgecolor="white",
        palette={"Nested LODO": "#A8C9C7", "Locked panel": "#34256F"},
        ax=ax_summary,
    )
    handles_summary, labels_summary = ax_summary.get_legend_handles_labels()
    ax_summary.legend(handles_summary[:2], labels_summary[:2], frameon=False, loc="upper left", ncol=2)
    for i, disease in enumerate(FIG1_DISEASES):
        nrow = disease_summary[disease_summary["disease"].astype(str) == disease].iloc[0]
        lrow = final_summary[final_summary["disease"].astype(str) == disease].iloc[0]
        ax_summary.plot([i - 0.32, i - 0.10], [nrow["median_auc"], nrow["median_auc"]], color=TEXT, lw=2.0, zorder=5)
        ax_summary.plot([i + 0.10, i + 0.32], [lrow["median_auc"], lrow["median_auc"]], color=TEXT, lw=2.0, zorder=5)
    ax_summary.axhline(0.5, color="#9CA3AF", lw=1.0, ls="--", zorder=0)
    ax_summary.set_ylim(0.45, 1.02)
    ax_summary.set_ylabel("AUC")
    ax_summary.set_xlabel("")
    ax_summary.set_title("Cohort-level AUC distributions across evaluation modes", loc="left", fontweight="bold", color=TEXT)
    ax_summary.grid(axis="y", color=GRID, lw=0.8)
    ax_summary.spines[["top", "right"]].set_visible(False)
    panel_label(ax_summary, "C")

    stability = (
        final_features.groupby("disease", as_index=False)
        .agg(
            final_k=("final_k", "first"),
            mean_selection_frequency=("nested_selection_frequency", "mean"),
            q25=("nested_selection_frequency", lambda x: float(np.percentile(x, 25))),
            q75=("nested_selection_frequency", lambda x: float(np.percentile(x, 75))),
            stable_features=("nested_selection_frequency", lambda x: int((x >= 0.5).sum())),
        )
    )
    stability["disease"] = pd.Categorical(stability["disease"], categories=FIG1_DISEASES, ordered=True)
    stability = stability.sort_values("disease")
    y_stab = np.arange(stability.shape[0])
    for yi, (_, row) in zip(y_stab, stability.iterrows()):
        disease = str(row["disease"])
        ax_stability.plot([row["q25"], row["q75"]], [yi, yi], color=FIG1_COLORS[disease], lw=3, alpha=0.80, solid_capstyle="round")
        ax_stability.scatter(row["mean_selection_frequency"], yi, s=78, color=FIG1_COLORS[disease], edgecolor="white", linewidth=0.8, zorder=3)
        ax_stability.text(
            min(float(row["mean_selection_frequency"]) + 0.04, 0.98),
            yi,
            f"k={int(row['final_k'])}; stable={int(row['stable_features'])}",
            va="center",
            ha="left",
            fontsize=8,
            color=TEXT,
        )
    ax_stability.set_yticks(y_stab)
    ax_stability.set_yticklabels(stability["disease"].astype(str))
    ax_stability.invert_yaxis()
    ax_stability.set_xlim(0, 1.05)
    ax_stability.set_xlabel("Nested selection frequency")
    ax_stability.set_title("Locked-panel stability summary", loc="left", fontweight="bold", color=TEXT)
    ax_stability.grid(axis="x", color=GRID, lw=0.8)
    ax_stability.spines[["top", "right"]].set_visible(False)
    panel_label(ax_stability, "E")

    y_final = np.arange(final_auc.shape[0])
    final_colors = [FIG1_COLORS[str(d)] for d in final_auc["disease"]]
    for yi, (_, row) in zip(y_final, final_auc.iterrows()):
        color = FIG1_COLORS[str(row["disease"])]
        ax_final_auc.plot([float(row["auc_ci_low"]), float(row["auc_ci_high"])], [yi, yi], color=color, lw=2.1, alpha=0.72, solid_capstyle="round", zorder=2)
        ax_final_auc.plot([float(row["auc_ci_low"]), float(row["auc_ci_low"])], [yi - 0.12, yi + 0.12], color=color, lw=1.3, alpha=0.72)
        ax_final_auc.plot([float(row["auc_ci_high"]), float(row["auc_ci_high"])], [yi - 0.12, yi + 0.12], color=color, lw=1.3, alpha=0.72)
    ax_final_auc.scatter(final_auc["auc"], y_final, s=58, c=final_colors, edgecolor="white", linewidth=0.8, zorder=3)
    ax_final_auc.axvline(0.5, color="#9CA3AF", lw=1.0, ls="--", zorder=0)
    for yi, (_, row) in zip(y_final, final_auc.iterrows()):
        sig = "*" if bool(row.get("significant_auc_gt_0_5", False)) else ""
        ax_final_auc.text(min(float(row["auc_ci_high"]) + 0.012, 1.055), yi, f"{float(row['auc']):.2f}{sig}", va="center", ha="left", fontsize=7.5, color=TEXT)
    for disease in FIG1_DISEASES:
        idx = final_auc.index[final_auc["disease"].astype(str) == disease].to_numpy()
        if idx.size:
            summary_row = final_summary[final_summary["disease"].astype(str) == disease].iloc[0]
            ax_final_auc.text(
                0.462,
                idx.min() - 0.48,
                f"{disease}: locked k={int(summary_row['final_effective_k'])}, median AUC={float(summary_row['median_auc']):.2f}",
                ha="left",
                va="center",
                fontsize=8.2,
                fontweight="bold",
                color=FIG1_COLORS[disease],
            )
    ax_final_auc.set_yticks(y_final)
    ax_final_auc.set_yticklabels(final_auc["row_label"])
    ax_final_auc.invert_yaxis()
    ax_final_auc.set_ylim(final_auc.shape[0] - 0.4, -1.05)
    ax_final_auc.set_xlim(0.45, 1.10)
    ax_final_auc.set_xlabel("AUC with bootstrap 95% CI")
    ax_final_auc.set_title("Locked signature transfer performance", loc="left", fontweight="bold", color=TEXT)
    ax_final_auc.grid(axis="x", color=GRID, lw=0.8)
    ax_final_auc.spines[["top", "right"]].set_visible(False)
    panel_label(ax_final_auc, "D")

    fig.legend(handles=handles, frameon=False, loc="upper center", ncol=3, bbox_to_anchor=(0.58, 0.995))
    fig.suptitle("Cross-cohort performance and locked species signatures", x=0.05, y=0.995, ha="left", fontsize=13, fontweight="bold", color=TEXT)
    fig.text(
        0.05,
        0.968,
        "A-C show strictly nested LODO validation; D-E show full-disease locked panels with target cohorts excluded from model training. Colors follow a light cyan-gray to deep indigo scale.",
        ha="left",
        va="top",
        fontsize=8.5,
        color=MUTED,
    )
    fig.savefig(OUTDIR / "figure1_lodo_bestk_auc.pdf", bbox_inches="tight")
    fig.savefig(OUTDIR / "figure1_lodo_bestk_auc.png", dpi=300, bbox_inches="tight")
    fig.savefig(OUTDIR / "figure1_integrated_refined.pdf", bbox_inches="tight")
    fig.savefig(OUTDIR / "figure1_integrated_refined.png", dpi=300, bbox_inches="tight")
    plt.close(fig)

    draw_figure1_top_species_supplement(final_features)


def draw_figure1_top_species_supplement(final_features: pd.DataFrame) -> None:
    plot_features = final_features.copy()
    plot_features["disease"] = pd.Categorical(plot_features["disease"], categories=FIG1_DISEASES, ordered=True)
    fig, axes = plt.subplots(1, 3, figsize=(15.5, 7.0), sharex=False)
    for ax, disease in zip(axes, FIG1_DISEASES):
        data = (
            plot_features[plot_features["disease"].astype(str) == disease]
            .sort_values("feature_rank")
            .head(15)
            .copy()
        )
        data = data.iloc[::-1].reset_index(drop=True)
        y = np.arange(data.shape[0])
        x = data["full_disease_log2fc"].astype(float).to_numpy()
        sizes = 38 + 150 * data["nested_selection_frequency"].astype(float).to_numpy()
        ax.axvline(0, color="#9CA3AF", lw=1.0, ls="--", zorder=0)
        for yi, xi in zip(y, x):
            ax.plot([0, xi], [yi, yi], color=FIG1_COLORS[disease], lw=1.4, alpha=0.45, zorder=1)
        ax.scatter(
            x,
            y,
            s=sizes,
            color=FIG1_COLORS[disease],
            edgecolor="white",
            linewidth=0.8,
            alpha=0.94,
            zorder=3,
        )
        labels = []
        for _, row in data.iterrows():
            name = clean_species(row["species"])
            if len(name) > 36:
                name = name[:33] + "..."
            labels.append(f"{int(row['feature_rank'])}. {name}")
        ax.set_yticks(y)
        ax.set_yticklabels(labels, fontsize=7)
        ax.set_title(f"{disease} locked panel top species", loc="left", fontweight="bold", color=FIG1_COLORS[disease])
        ax.set_xlabel("Full-disease log2FC")
        ax.grid(axis="x", color=GRID, lw=0.8)
        ax.spines[["top", "right"]].set_visible(False)
    size_handles = [
        plt.scatter([], [], s=38 + 150 * value, color="#6F7DBA", edgecolor="white", linewidth=0.8, label=f"{value:.2f}")
        for value in (0.25, 0.50, 0.75, 1.00)
    ]
    fig.legend(handles=size_handles, title="Nested selection frequency", frameon=False, loc="upper center", ncol=4, bbox_to_anchor=(0.58, 0.98))
    fig.suptitle(
        "Supplementary Figure 1. Top species in final locked disease panels",
        x=0.04,
        y=0.995,
        ha="left",
        fontsize=13,
        fontweight="bold",
        color=TEXT,
    )
    fig.text(
        0.04,
        0.955,
        "Lollipop position shows full-disease log2 fold change; point size shows how often each species was selected across nested LODO folds.",
        ha="left",
        va="top",
        fontsize=8.5,
        color=MUTED,
    )
    fig.savefig(OUTDIR / "figure1_supplement_top_species.pdf", bbox_inches="tight")
    fig.savefig(OUTDIR / "figure1_supplement_top_species.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def membership_counts(assoc: pd.DataFrame) -> pd.DataFrame:
    memberships = []
    for species, group in assoc.groupby("species"):
        diseases = tuple(d for d in DISEASES if d in set(group["disease"]))
        memberships.append({"species": species, "diseases": "|".join(diseases), "n_diseases": len(diseases)})
    df = pd.DataFrame(memberships)
    counts = df.groupby(["diseases", "n_diseases"], as_index=False).size().rename(columns={"size": "species_count"})
    return counts.sort_values(["n_diseases", "species_count", "diseases"], ascending=[False, False, True]).reset_index(drop=True)


def draw_upset(ax_bar: plt.Axes, ax_matrix: plt.Axes, counts: pd.DataFrame) -> None:
    counts = counts.copy()
    counts = counts[counts["n_diseases"] >= 2].head(12).reset_index(drop=True)
    x = np.arange(counts.shape[0])
    bar_colors = ["#111827" if n == 4 else "#475569" if n == 3 else "#94A3B8" for n in counts["n_diseases"]]
    ax_bar.bar(x, counts["species_count"], color=bar_colors, width=0.72)
    for xi, count in zip(x, counts["species_count"]):
        ax_bar.text(xi, count + 2, str(int(count)), ha="center", va="bottom", fontsize=8, color=TEXT)
    ax_bar.set_ylabel("Species")
    ax_bar.set_title("Top shared intersections among disease-associated species", loc="left", fontweight="bold", color=TEXT)
    ax_bar.grid(axis="y", color=GRID, lw=0.8)
    ax_bar.spines[["top", "right"]].set_visible(False)
    ax_bar.set_xticks([])
    panel_label(ax_bar, "A")

    ax_matrix.set_xlim(-0.5, counts.shape[0] - 0.5)
    ax_matrix.set_ylim(-0.5, len(DISEASES) - 0.5)
    for xi, disease_set in enumerate(counts["diseases"]):
        members = set(str(disease_set).split("|"))
        selected_y = []
        for yi, disease in enumerate(DISEASES[::-1]):
            y = yi
            selected = disease in members
            ax_matrix.scatter(
                xi,
                y,
                s=58 if selected else 24,
                color=DISEASE_COLORS[disease] if selected else "#E5E7EB",
                edgecolor="white",
                linewidth=0.6,
                zorder=3,
            )
            if selected:
                selected_y.append(y)
        if len(selected_y) > 1:
            ax_matrix.plot([xi, xi], [min(selected_y), max(selected_y)], color="#334155", lw=1.0, zorder=2)
    ax_matrix.set_yticks(np.arange(len(DISEASES)))
    ax_matrix.set_yticklabels(DISEASES[::-1])
    ax_matrix.set_xticks(x)
    ax_matrix.set_xticklabels(counts["diseases"].str.replace("|", "+", regex=False), rotation=45, ha="right")
    ax_matrix.tick_params(axis="x", length=0)
    ax_matrix.tick_params(axis="y", length=0)
    for spine in ax_matrix.spines.values():
        spine.set_visible(False)


def draw_pairwise_heatmap(ax: plt.Axes, overlap: pd.DataFrame) -> None:
    pair = overlap[overlap["n_diseases"] == 2].copy()
    ratio = pd.DataFrame(np.nan, index=DISEASES, columns=DISEASES, dtype=float)
    annot = pd.DataFrame("", index=DISEASES, columns=DISEASES)
    for _, row in pair.iterrows():
        d1, d2 = str(row["diseases"]).split("|")
        for a, b in [(d1, d2), (d2, d1)]:
            ratio.loc[a, b] = float(row["observed_expected_ratio"])
            annot.loc[a, b] = f"{int(row['shared_count'])}\n{float(row['observed_expected_ratio']):.2f}x\n{q_label(float(row['permutation_q']))}"
    np.fill_diagonal(ratio.values, 1.0)
    for disease in DISEASES:
        annot.loc[disease, disease] = "-"
    sns.heatmap(
        ratio,
        ax=ax,
        cmap="YlGnBu",
        vmin=0.8,
        vmax=max(2.0, np.nanmax(ratio.to_numpy())),
        annot=annot,
        fmt="",
        linewidths=0.6,
        linecolor="white",
        cbar_kws={"label": "Observed / expected", "shrink": 0.75},
    )
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_title("Pairwise overlap is enriched above a matched detectable-species background", loc="left", fontweight="bold", color=TEXT)
    panel_label(ax, "B")


def draw_multiway(ax: plt.Axes, overlap: pd.DataFrame) -> None:
    multi = overlap[overlap["n_diseases"] >= 3].copy()
    multi = multi.sort_values(["n_diseases", "shared_count"], ascending=[True, False]).reset_index(drop=True)
    x = np.arange(multi.shape[0])
    colors = ["#111827" if n == 4 else "#64748B" for n in multi["n_diseases"]]
    ax.bar(x, multi["shared_count"], color=colors, width=0.72, label="Observed")
    ax.scatter(x, multi["expected_permutation"], color="#F59E0B", s=30, zorder=4, label="Expected")
    for xi, (_, row) in zip(x, multi.iterrows()):
        ax.text(xi, float(row["shared_count"]) + 0.8, q_label(float(row["permutation_q"])), ha="center", va="bottom", fontsize=7, color=TEXT)
    ax.set_xticks(x)
    ax.set_xticklabels(multi["diseases"].str.replace("|", "\n", regex=False))
    ax.set_ylabel("Shared species")
    ax.set_title("Multi-disease intersections remain non-random", loc="left", fontweight="bold", color=TEXT)
    ax.grid(axis="y", color=GRID, lw=0.8)
    ax.legend(frameon=False, loc="upper right")
    ax.spines[["top", "right"]].set_visible(False)
    panel_label(ax, "C")


def draw_shared_heatmap(ax: plt.Axes, shared: pd.DataFrame) -> None:
    plot = shared[(shared["n_diseases"] >= 3) & (shared["named_species"] == "yes")].copy()
    score_cols = [f"{d}_rank_score" for d in DISEASES]
    plot["total_rank_score"] = plot[score_cols].fillna(0).sum(axis=1)
    plot = plot.sort_values(["n_diseases", "same_direction", "total_rank_score"], ascending=[False, False, False]).head(18)
    effect_cols = [f"{d}_meta_log2fc" for d in DISEASES]
    matrix = plot.set_index("species")[effect_cols]
    matrix.columns = DISEASES
    matrix.index = [clean_species(i) for i in matrix.index]
    sns.heatmap(
        matrix,
        ax=ax,
        cmap="vlag",
        center=0,
        vmin=-3,
        vmax=3,
        linewidths=0.4,
        linecolor="white",
        cbar_kws={"label": "Meta log2FC", "shrink": 0.70},
    )
    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_title("Representative shared species preserve disease-specific effect directions", loc="left", fontweight="bold", color=TEXT)
    ax.tick_params(axis="y", labelsize=7)
    ax.tick_params(axis="x", rotation=0)
    panel_label(ax, "D")


def write_figure2_source_data(assoc: pd.DataFrame, overlap: pd.DataFrame, counts: pd.DataFrame) -> None:
    membership = counts.copy()
    membership.insert(0, "panel", "upset_membership")
    overlap_out = overlap.copy()
    overlap_out.insert(0, "panel", "overlap_enrichment")
    common_cols = sorted(set(membership.columns).union(overlap_out.columns))
    out = pd.concat(
        [membership.reindex(columns=common_cols), overlap_out.reindex(columns=common_cols)],
        ignore_index=True,
    )
    out.to_csv(OUTDIR / "figure2_source_data.tsv", sep="\t", index=False)


def draw_figure2() -> None:
    assoc = pd.read_csv(ASSOC_DIR / "disease_associated_species_top300.tsv", sep="\t")
    overlap = pd.read_csv(ASSOC_DIR / "overlap_summary_top300.tsv", sep="\t")
    shared = pd.read_csv(ASSOC_DIR / "shared_associated_species_top300.tsv", sep="\t")
    counts = membership_counts(assoc)
    write_figure2_source_data(assoc, overlap, counts)

    fig = plt.figure(figsize=(13.2, 10.2))
    gs = gridspec.GridSpec(3, 2, height_ratios=[1.35, 0.58, 2.2], hspace=0.38, wspace=0.30, figure=fig)
    ax_upset_bar = fig.add_subplot(gs[0, 0])
    ax_upset_matrix = fig.add_subplot(gs[1, 0])
    ax_pairwise = fig.add_subplot(gs[0:2, 1])
    ax_multiway = fig.add_subplot(gs[2, 0])
    ax_heatmap = fig.add_subplot(gs[2, 1])

    draw_upset(ax_upset_bar, ax_upset_matrix, counts)
    draw_pairwise_heatmap(ax_pairwise, overlap)
    draw_multiway(ax_multiway, overlap)
    draw_shared_heatmap(ax_heatmap, shared)

    set_sizes = assoc.groupby("disease")["species"].nunique().reindex(DISEASES).to_dict()
    shared2 = assoc.groupby("species")["disease"].nunique()
    n_shared2 = int((shared2 >= 2).sum())
    n_shared3 = int((shared2 >= 3).sum())
    n_shared4 = int((shared2 == 4).sum())
    subtitle = (
        "Top-ranked disease-associated species sets "
        + ", ".join(f"{d}=Top {int(set_sizes[d])}" for d in DISEASES)
        + f"; shared by >=2 diseases={n_shared2}, >=3 diseases={n_shared3}, all four={n_shared4}."
    )
    fig.suptitle(
        "Disease-associated gut species show extensive overlap across CRA, CRC, IBD and IBS",
        x=0.04,
        y=0.992,
        ha="left",
        fontsize=13,
        fontweight="bold",
        color=TEXT,
    )
    fig.text(0.04, 0.963, subtitle, ha="left", va="top", fontsize=8.5, color=MUTED)
    fig.savefig(OUTDIR / "figure2_associated_species_overlap_top300.pdf", bbox_inches="tight")
    fig.savefig(OUTDIR / "figure2_associated_species_overlap_top300.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def write_legends() -> None:
    lines = [
        "# Figure legends",
        "",
        "## Figure 1. Cross-cohort performance and locked species signatures",
        "",
        "Leave-one-dataset-out validation was performed independently for CRA, CRC and IBD. "
        "For each held-out dataset, all prevalence and abundance filtering, disease-association ranking, "
        "correlation-based feature de-redundancy, best-k selection and standardization were estimated only from the training datasets. "
        "A balanced logistic regression classifier was then evaluated in the unseen cohort. "
        "Panel A shows held-out AUCs with stratified bootstrap 95% confidence intervals and annotates the training-selected best-k feature count for each fold. "
        "Panel B shows the disease-level best-k tuning landscape. Panel C compares nested LODO and locked-panel cohort-level AUC distributions. "
        "Panel D shows locked-panel transfer AUCs, with final k and disease-level median AUC annotated directly on the panel. "
        "Panel E summarizes locked-panel stability using nested feature-selection frequencies. Locked-panel results are post-selection transfer estimates, whereas Panel A remains the strictly nested no-test-leakage validation.",
        "",
        "## Supplementary Figure 1. Top species in final locked disease panels",
        "",
        "Top species are visualized as lollipop/dot plots. The x-axis reports the full-disease log2 fold change, and point size reports how often each species was selected across nested LODO folds.",
        "",
        "## Figure 2. Disease-associated species overlap across CRA, CRC, IBD and IBS",
        "",
        "The overlap analysis used the Top300 ranked disease-associated species for each disease, not the stricter Top50 key-species set. "
        "Panel A shows the leading shared intersections, Panel B reports pairwise observed/expected enrichment against a matched detectable-species background, "
        "Panel C compares observed and permutation-expected multi-disease intersections, and Panel D shows representative shared species with disease-level meta log2 fold changes. "
        "This figure supports the conclusion that disease-associated bacteria display extensive species-level overlap across CRA, CRC, IBD and IBS.",
        "",
    ]
    (OUTDIR / "figure_legends.md").write_text("\n".join(lines), encoding="utf-8")


def validate_outputs() -> None:
    expected = [
        "figure1_lodo_bestk_auc.pdf",
        "figure1_lodo_bestk_auc.png",
        "figure1_integrated_refined.pdf",
        "figure1_integrated_refined.png",
        "figure1_supplement_top_species.pdf",
        "figure1_supplement_top_species.png",
        "figure2_associated_species_overlap_top300.pdf",
        "figure2_associated_species_overlap_top300.png",
        "figure1_source_data.tsv",
        "figure1_selected_features.tsv",
        "figure1_predictions.tsv",
        "figure1_disease_summary.tsv",
        "figure1_final_signature_features.tsv",
        "figure1_final_signature_auc_by_cohort.tsv",
        "figure1_final_signature_predictions.tsv",
        "figure1_final_signature_disease_summary.tsv",
        "figure2_source_data.tsv",
        "figure_legends.md",
    ]
    missing = []
    for name in expected:
        path = OUTDIR / name
        if not path.exists() or path.stat().st_size == 0:
            missing.append(name)
    if missing:
        raise RuntimeError(f"Missing or empty outputs: {missing}")


def run() -> None:
    OUTDIR.mkdir(parents=True, exist_ok=True)
    setup_style()
    draw_figure1()
    draw_figure2()
    write_legends()
    validate_outputs()
    print(f"Wrote requested figures to {OUTDIR}", flush=True)


if __name__ == "__main__":
    run()
