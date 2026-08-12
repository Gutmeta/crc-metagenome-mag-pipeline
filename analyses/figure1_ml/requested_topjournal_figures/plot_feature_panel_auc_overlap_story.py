#!/usr/bin/env python3
"""Three-panel publication figure for locked panels and IBS-related overlap."""

from __future__ import annotations

import os
from collections import Counter
from itertools import combinations
from pathlib import Path

os.environ.setdefault("MPLCONFIGDIR", str(Path(__file__).resolve().parent / ".mplconfig"))

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
import pandas as pd
from matplotlib.gridspec import GridSpec, GridSpecFromSubplotSpec
from scipy.stats import gaussian_kde


ROOT = Path(__file__).resolve().parents[1]
OUTDIR = ROOT / "requested_topjournal_figures"
ASSOC_DIR = ROOT / "species_association_overlap"

FIG_DISEASES = ("CRA", "CRC", "IBD")
OVERLAP_SETS = ("CRA", "CRC", "IBD", "IBS")
COLORS = {
    "CRA": "#0B3C68",
    "CRC": "#2F7EA8",
    "IBD": "#8ABBD3",
    "IBS": "#E9BFAE",
}
DARK_COLORS = {
    "CRA": "#08365F",
    "CRC": "#1F668D",
    "IBD": "#5B9BBC",
    "IBS": "#C98773",
}
AUC_CMAP = LinearSegmentedColormap.from_list(
    "auc_reference_blue_peach",
    ["#F1C7B5", "#E9F0F4", "#8EC1DC", "#2F7EA8", "#0B3C68"],
)
PAIRWISE_CMAP = LinearSegmentedColormap.from_list(
    "reference_blue_enrichment",
    ["#F7FAFC", "#E9F0F4", "#BFD9E6", "#6BA7C6", "#0B3C68"],
)
TEXT = "#111820"
MUTED = "#52616E"
GRID = "#D9E1E8"
BAND = "#F7FAFC"
ZERO_LINE = "#8B98A5"
SPINE = "#1F2933"
EXPECTED = "#D7A187"
CONTROL_FILL = "#7897AC"
CASE_FILL = "#E4A38E"
CONTROL_EDGE = "#274C63"
CASE_EDGE = "#9D5B4B"
SIG_BAR = "#0B3C68"
NONSIG_BAR = "#BFD9E6"
FIGURE_LEFT_ANCHOR = 0.09126
PANEL_A_TOP = 0.935


def setup_style() -> None:
    plt.rcParams.update(
        {
            "font.family": "Arial",
            "font.sans-serif": ["Arial"],
            "font.size": 11.4,
            "font.weight": "normal",
            "axes.titlesize": 12.4,
            "axes.titleweight": "normal",
            "axes.labelsize": 10.8,
            "axes.labelweight": "normal",
            "xtick.labelsize": 9.8,
            "ytick.labelsize": 9.4,
            "axes.edgecolor": "#CBD5E1",
            "axes.linewidth": 0.8,
            "axes.labelcolor": TEXT,
            "xtick.color": TEXT,
            "ytick.color": TEXT,
            "figure.facecolor": "white",
            "axes.facecolor": "white",
            "savefig.facecolor": "white",
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
        }
    )


def panel_label(ax: plt.Axes, label: str, x: float = -0.12, y: float = 1.04) -> None:
    ax.text(
        x,
        y,
        label,
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=17.0,
        fontweight="normal",
        color=TEXT,
    )


def truncate_label(text: str, max_len: int = 25) -> str:
    text = str(text)
    return text if len(text) <= max_len else text[: max_len - 1] + "..."


COHORT_LABEL_OVERRIDES = {
    "Damman_2015_FMT_UC": "Damman 2015",
    "Bushman_2020_pediatric_Cdiff_IBD": "Bushman 2020",
    "Lloyd-Price_2019_HMP2IBD": "Lloyd-Price 2019",
    "Douglas_2018_child_Crohn": "Douglas 2018",
}


def clean_cohort(name: str) -> str:
    text = str(name)
    if text in COHORT_LABEL_OVERRIDES:
        return COHORT_LABEL_OVERRIDES[text]
    remove_tokens = [
        "_CRC_",
        "_IBD_",
        "_UC_",
        "_Crohn",
        "_CRC",
        "_IBD",
        "_China",
        "_Italy",
        "_France",
        "_Austria",
        "_Germany",
        "_USA",
        "_pediatric",
    ]
    for token in remove_tokens:
        text = text.replace(token, "_")
    text = text.strip("_").replace("_", " ")
    return " ".join(text.split())


def q_text(q: float) -> str:
    if not np.isfinite(q):
        return "n.a."
    if q < 0.001:
        return "q<0.001"
    if q < 0.01:
        return f"q={q:.3f}"
    return f"q={q:.2f}"


def load_inputs() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    auc = pd.read_csv(OUTDIR / "figure1_final_signature_auc_by_cohort.tsv", sep="\t")
    auc = auc[auc["status"].eq("ok")].copy()
    auc["disease"] = pd.Categorical(auc["disease"], categories=FIG_DISEASES, ordered=True)
    auc = auc.sort_values(["disease", "auc", "left_out_cohort"], ascending=[True, True, True]).reset_index(drop=True)
    if auc.empty:
        raise RuntimeError("No successful locked-panel AUC rows were found.")

    predictions = pd.read_csv(OUTDIR / "figure1_final_signature_predictions.tsv", sep="\t")
    predictions = predictions[predictions["disease"].isin(FIG_DISEASES)].copy()
    predictions["predicted_probability"] = pd.to_numeric(predictions["predicted_probability"], errors="coerce")
    predictions = predictions.dropna(subset=["predicted_probability"])

    final_features = pd.read_csv(OUTDIR / "figure1_final_signature_features.tsv", sep="\t")
    assoc = pd.read_csv(ASSOC_DIR / "disease_associated_species_top300.tsv", sep="\t")
    overlap = pd.read_csv(ASSOC_DIR / "overlap_summary_top300.tsv", sep="\t")
    return auc, predictions, final_features, assoc, overlap


def associated_sets(assoc: pd.DataFrame) -> dict[str, set[str]]:
    return {
        disease: set(assoc.loc[assoc["disease"].eq(disease), "species"].astype(str))
        for disease in OVERLAP_SETS
    }


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


def auc_color(value: float) -> tuple[float, float, float, float]:
    normed = float(np.clip((value - 0.48) / 0.36, 0, 1))
    return AUC_CMAP(normed)


def density_curve(values: np.ndarray, grid: np.ndarray) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    values = values[np.isfinite(values)]
    if values.size < 3 or np.nanstd(values) < 1e-4:
        center = float(np.nanmean(values)) if values.size else 0.5
        dens = np.exp(-0.5 * ((grid - center) / 0.035) ** 2)
    else:
        try:
            dens = gaussian_kde(values, bw_method=0.28)(grid)
        except Exception:
            center = float(np.nanmean(values))
            dens = np.exp(-0.5 * ((grid - center) / 0.05) ** 2)
    max_density = float(np.nanmax(dens)) if np.size(dens) else 0.0
    return dens / max_density if max_density > 0 else dens


def draw_auc_forest(
    fig: plt.Figure,
    spec: GridSpecFromSubplotSpec,
    auc: pd.DataFrame,
    predictions: pd.DataFrame,
) -> pd.DataFrame:
    panel = GridSpecFromSubplotSpec(1, 3, subplot_spec=spec, wspace=0.28)
    rows: list[dict[str, object]] = []

    for idx, disease in enumerate(FIG_DISEASES):
        subpanel = GridSpecFromSubplotSpec(
            1,
            2,
            subplot_spec=panel[0, idx],
            width_ratios=[0.84, 1.24],
            wspace=0.09,
        )
        ax_auc = fig.add_subplot(subpanel[0, 0])
        ax_prob = fig.add_subplot(subpanel[0, 1], sharey=ax_auc)
        sub_auc = auc[auc["disease"].astype(str).eq(disease)].sort_values("auc", ascending=False).reset_index(drop=True)
        y = np.arange(sub_auc.shape[0], dtype=float)
        dark = DARK_COLORS[disease]

        for yy in y:
            if int(yy) % 2 == 0:
                ax_auc.axhspan(yy - 0.48, yy + 0.48, color=BAND, zorder=0)
                ax_prob.axhspan(yy - 0.48, yy + 0.48, color=BAND, zorder=0)
        ax_auc.axvline(0.5, color=ZERO_LINE, lw=0.9, ls=(0, (4, 3)), zorder=1)
        ax_prob.axvline(0.5, color=ZERO_LINE, lw=0.9, ls=(0, (4, 3)), zorder=1)

        for yy, (_, row) in zip(y, sub_auc.iterrows()):
            cohort = str(row["left_out_cohort"])
            pred = predictions[
                predictions["disease"].astype(str).eq(disease)
                & predictions["left_out_cohort"].astype(str).eq(cohort)
            ]
            control_pred = pred.loc[pred["y_true"].eq(0), "predicted_probability"].to_numpy(dtype=float)
            case_pred = pred.loc[pred["y_true"].eq(1), "predicted_probability"].to_numpy(dtype=float)
            auc_value = float(row["auc"])
            ci_low = float(row["auc_ci_low"])
            ci_high = float(row["auc_ci_high"])
            fill = auc_color(auc_value)

            ci_low_plot = max(ci_low, 0.45)
            ci_high_plot = min(ci_high, 1.0)
            xerr = np.array(
                [
                    [max(auc_value - ci_low_plot, 0.0)],
                    [max(ci_high_plot - auc_value, 0.0)],
                ]
            )
            ax_auc.errorbar(
                auc_value,
                yy,
                xerr=xerr,
                fmt="o",
                markersize=4.8,
                markerfacecolor=fill,
                markeredgecolor=dark,
                markeredgewidth=0.62,
                ecolor="#83919E",
                elinewidth=1.25,
                capsize=2.4,
                capthick=0.95,
                zorder=3,
            )

            control_q25, control_med, control_q75 = np.nanpercentile(control_pred, [25, 50, 75])
            case_q25, case_med, case_q75 = np.nanpercentile(case_pred, [25, 50, 75])
            y_case = yy - 0.16
            y_control = yy + 0.16
            ax_prob.plot(
                [control_q25, control_q75],
                [y_control, y_control],
                color=CONTROL_EDGE,
                lw=3.1,
                alpha=0.96,
                solid_capstyle="round",
                zorder=2,
            )
            ax_prob.scatter(control_med, y_control, s=27, color=CONTROL_FILL, edgecolor=CONTROL_EDGE, linewidth=0.55, zorder=3)
            ax_prob.plot(
                [case_q25, case_q75],
                [y_case, y_case],
                color=CASE_EDGE,
                lw=3.1,
                alpha=0.96,
                solid_capstyle="round",
                zorder=2,
            )
            ax_prob.scatter(case_med, y_case, s=27, color=CASE_FILL, edgecolor=CASE_EDGE, linewidth=0.55, zorder=3)
            ax_prob.plot(
                [control_med, case_med],
                [yy, yy],
                color="#CBD5DF",
                lw=0.55,
                alpha=0.85,
                zorder=1,
            )
            rows.append(
                {
                    **row.to_dict(),
                    "plot_y": float(yy),
                    "label": clean_cohort(cohort),
                    "control_n": int(control_pred.size),
                    "case_n": int(case_pred.size),
                    "control_median_probability": float(control_med),
                    "control_q25_probability": float(control_q25),
                    "control_q75_probability": float(control_q75),
                    "case_median_probability": float(case_med),
                    "case_q25_probability": float(case_q25),
                    "case_q75_probability": float(case_q75),
                }
            )

        labels = [
            f"{truncate_label(clean_cohort(row['left_out_cohort']), 16)}\n({int(row['test_case'])}/{int(row['test_control'])})"
            for _, row in sub_auc.iterrows()
        ]
        final_k = int(sub_auc["final_k"].iloc[0])
        median_auc = float(sub_auc["auc"].median())
        ax_auc.set_title(f"{disease} species", loc="left", color=TEXT, fontweight="normal", pad=7, fontsize=12.2)
        ax_prob.text(
            0.98,
            1.012,
            f"best K={final_k} | median AUC={median_auc:.2f}",
            transform=ax_prob.transAxes,
            ha="right",
            va="bottom",
            fontsize=8.6,
            color=dark,
            fontweight="normal",
            bbox={
                "boxstyle": "round,pad=0.14,rounding_size=0.06",
                "facecolor": "#F2F7FA",
                "edgecolor": "#D9E7EF",
                "linewidth": 0.55,
            },
        )
        ax_auc.set_xlim(0.45, 1.0)
        ax_auc.set_xticks([0.5, 0.75, 1.0])
        ax_auc.set_xticklabels(["0.5", "0.75", "1"])
        ax_auc.set_xlabel("AUC", labelpad=4)
        ax_auc.set_ylim(sub_auc.shape[0] - 0.62, -0.86)
        ax_auc.set_yticks(y)
        ax_auc.set_yticklabels(labels, fontsize=8.5, linespacing=0.96)
        ax_prob.set_xlim(0.0, 1.0)
        ax_prob.set_xticks([0.0, 0.5, 1.0])
        ax_prob.set_xticklabels(["0", "0.5", "1"])
        ax_prob.set_xlabel("Predicted probability", labelpad=4)
        ax_prob.tick_params(axis="y", labelleft=False, length=0)
        if idx == 0:
            fig.text(
                0.035,
                PANEL_A_TOP,
                "A",
                ha="left",
                va="top",
                fontsize=19.2,
                fontweight="normal",
                color=TEXT,
            )
            handles = [
                plt.Line2D([0], [0], color=CONTROL_EDGE, lw=2.4, marker="o", markersize=3.5, label="Control IQR/median"),
                plt.Line2D([0], [0], color=CASE_EDGE, lw=2.4, marker="o", markersize=3.5, label="Case IQR/median"),
            ]
            fig.legend(
                handles=handles,
                frameon=False,
                loc="lower left",
                bbox_to_anchor=(FIGURE_LEFT_ANCHOR - 0.0024, 0.012),
                ncol=2,
                fontsize=8.8,
                handlelength=1.55,
                borderpad=0.1,
                labelspacing=0.3,
                columnspacing=1.4,
            )
        ax_auc.tick_params(axis="y", length=0, pad=2)
        for axis in (ax_auc, ax_prob):
            axis.set_axisbelow(True)
            axis.grid(axis="x", color=GRID, lw=0.58, alpha=0.78)
            axis.tick_params(axis="x", length=3, width=0.8, color=SPINE)
            for spine in axis.spines.values():
                spine.set_linewidth(0.85)
                spine.set_color(SPINE)
            axis.spines["top"].set_visible(False)
            axis.spines["right"].set_visible(False)
        ax_prob.spines["left"].set_visible(False)

    return pd.DataFrame(rows)


def draw_upset(ax_bar: plt.Axes, ax_matrix: plt.Axes, counts: pd.DataFrame) -> None:
    counts = counts[counts["n_diseases"].ge(2)].head(10).reset_index(drop=True)
    x = np.arange(counts.shape[0])
    bar_colors = ["#111827" if n == 4 else "#475569" if n == 3 else "#94A3B8" for n in counts["n_diseases"]]
    ax_bar.bar(x, counts["species_count"], color=bar_colors, width=0.68)
    for xi, count in zip(x, counts["species_count"]):
        ax_bar.text(xi, count + 1.5, str(int(count)), ha="center", va="bottom", fontsize=10.1, color=TEXT)
    ax_bar.set_ylabel("No. species")
    ax_bar.set_ylim(0, max(counts["species_count"]) * 1.25)
    ax_bar.set_title("Disease-associated feature sharing", loc="left", fontweight="normal", color=TEXT, pad=8)
    ax_bar.grid(axis="y", color=GRID, lw=0.75)
    ax_bar.spines[["top", "right"]].set_visible(False)
    ax_bar.set_xticks([])
    panel_label(ax_bar, "B", x=-0.09, y=1.12)

    ax_matrix.set_xlim(-0.5, counts.shape[0] - 0.5)
    ax_matrix.set_ylim(-0.5, len(OVERLAP_SETS) - 0.5)
    for xi, combo in enumerate(counts["diseases"]):
        members = set(str(combo).split("|"))
        selected_y = []
        for yi, disease in enumerate(OVERLAP_SETS[::-1]):
            selected = disease in members
            ax_matrix.scatter(
                xi,
                yi,
                s=48 if selected else 20,
                color=COLORS[disease] if selected else "#E9EDF2",
                edgecolor="white",
                linewidth=0.55,
                zorder=3,
            )
            if selected:
                selected_y.append(yi)
        if len(selected_y) > 1:
            ax_matrix.plot([xi, xi], [min(selected_y), max(selected_y)], color="#334155", lw=0.9, zorder=2)
    ax_matrix.set_yticks(np.arange(len(OVERLAP_SETS)))
    ax_matrix.set_yticklabels(OVERLAP_SETS[::-1])
    ax_matrix.set_xticks(x)
    ax_matrix.set_xticklabels(counts["diseases"].str.replace("|", "+", regex=False), rotation=43, ha="right", fontsize=9.6)
    ax_matrix.tick_params(axis="both", length=0)
    for spine in ax_matrix.spines.values():
        spine.set_visible(False)


def draw_pairwise_heatmap(ax: plt.Axes, overlap: pd.DataFrame, sets: dict[str, set[str]]) -> pd.DataFrame:
    pair_rows = []
    for a, b in combinations(OVERLAP_SETS, 2):
        row = overlap[overlap["diseases"].eq(f"{a}|{b}")]
        if row.empty:
            row = overlap[overlap["diseases"].eq(f"{b}|{a}")]
        if row.empty:
            continue
        pair_rows.append(row.iloc[0].to_dict())

    ratio = pd.DataFrame(np.nan, index=OVERLAP_SETS, columns=OVERLAP_SETS)
    annot = pd.DataFrame("", index=OVERLAP_SETS, columns=OVERLAP_SETS)
    for row in pair_rows:
        a, b = str(row["diseases"]).split("|")
        shared = int(row["shared_count"])
        observed_expected = float(row["observed_expected_ratio"])
        for x, y in ((a, b), (b, a)):
            ratio.loc[x, y] = observed_expected
            annot.loc[x, y] = f"{shared}\n{observed_expected:.2g}x"
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
            weight = "normal"
            ax.text(
                j,
                i,
                annot.loc[disease_i, disease_j],
                ha="center",
                va="center",
                fontsize=9.0,
                color=color,
                fontweight=weight,
                linespacing=1.04,
            )
    ax.set_title("Pairwise enrichment", loc="left", fontweight="normal", color=TEXT, pad=6)
    ax.set_box_aspect(0.68)
    ax.set_anchor("W")
    for spine in ax.spines.values():
        spine.set_visible(False)
    panel_label(ax, "B", x=-0.18, y=1.08)
    return pd.DataFrame(pair_rows)


def draw_multiway(ax: plt.Axes, overlap: pd.DataFrame) -> pd.DataFrame:
    multi = overlap[(overlap["n_diseases"].ge(3)) & (overlap["n_diseases"].lt(4))].copy()
    combo_order = {combo: i for i, combo in enumerate(["CRA|CRC|IBD", "CRA|CRC|IBS", "CRA|IBD|IBS", "CRC|IBD|IBS"])}
    multi["combo_order"] = multi["diseases"].map(combo_order).fillna(99)
    multi = multi.sort_values(["combo_order", "shared_count"], ascending=[True, False]).reset_index(drop=True)
    x = np.arange(multi.shape[0])
    colors = [SIG_BAR if bool(sig) else NONSIG_BAR for sig in multi["significantly_above_random"]]
    ax.bar(x, multi["shared_count"], color=colors, edgecolor="white", linewidth=0.85, width=0.68, label="Observed", zorder=3)
    ax.scatter(x, multi["expected_permutation"], color=EXPECTED, s=42, edgecolor="white", linewidth=0.55, zorder=4, label="Expected")
    for xi, (_, row) in zip(x, multi.iterrows()):
        label = f"{float(row['observed_expected_ratio']):.1f}x"
        ax.text(xi, float(row["shared_count"]) + 0.75, label, ha="center", va="bottom", fontsize=8.8, color=TEXT)
    ax.set_xticks(x)
    xticklabels = ["+".join(str(combo).split("|")[:2]) + "\n+" + str(combo).split("|")[2] for combo in multi["diseases"]]
    ax.set_xticklabels(xticklabels, rotation=0, ha="center", fontsize=7.7, linespacing=0.95)
    ax.set_ylim(0, max(multi["shared_count"]) * 1.30)
    ax.set_ylabel("Shared species", labelpad=4, fontsize=10.0)
    ax.set_title("Multi-disease enrichment", loc="left", fontweight="normal", color=TEXT, pad=6)
    ax.set_axisbelow(True)
    ax.grid(axis="y", color=GRID, lw=0.75, alpha=0.9)
    ax.legend(frameon=False, loc="upper right", bbox_to_anchor=(1.0, 1.03), ncol=2, handlelength=1.05, borderpad=0.1, columnspacing=0.8, fontsize=8.8)
    ax.set_box_aspect(0.68)
    ax.set_anchor("W")
    panel_label(ax, "C", x=-0.22, y=1.16)
    ax.spines[["top", "right"]].set_visible(False)
    ax.spines["left"].set_color(SPINE)
    ax.spines["bottom"].set_color("#CBD5E1")
    return multi.drop(columns=["combo_order"], errors="ignore")


def draw_overlap_panel(fig: plt.Figure, spec: GridSpecFromSubplotSpec, counts: pd.DataFrame, overlap: pd.DataFrame, sets: dict[str, set[str]]) -> tuple[pd.DataFrame, pd.DataFrame]:
    panel = GridSpecFromSubplotSpec(2, 1, subplot_spec=spec, height_ratios=[1.0, 1.0], hspace=0.30)

    ax_heatmap = fig.add_subplot(panel[0])
    ax_multi = fig.add_subplot(panel[1])

    pairwise = draw_pairwise_heatmap(ax_heatmap, overlap, sets)
    multi = draw_multiway(ax_multi, overlap)
    return pairwise, multi


def locked_panel_overlap_audit(final_features: pd.DataFrame, assoc: pd.DataFrame) -> pd.DataFrame:
    sets: dict[str, set[str]] = {
        disease: set(final_features.loc[final_features["disease"].eq(disease), "species"].astype(str))
        for disease in FIG_DISEASES
    }
    sets["IBS-related"] = set(assoc.loc[assoc["disease"].eq("IBS"), "species"].astype(str))
    rows = []
    for a, b in combinations(sets, 2):
        shared = sets[a] & sets[b]
        rows.append(
            {
                "set_a": a,
                "set_b": b,
                "set_a_size": len(sets[a]),
                "set_b_size": len(sets[b]),
                "shared_count": len(shared),
                "overlap_coefficient": len(shared) / min(len(sets[a]), len(sets[b])) if min(len(sets[a]), len(sets[b])) else np.nan,
                "jaccard": len(shared) / len(sets[a] | sets[b]) if (sets[a] | sets[b]) else np.nan,
                "shared_species": ";".join(sorted(shared)),
            }
        )
    return pd.DataFrame(rows)


def prediction_distribution_audit(predictions: pd.DataFrame) -> pd.DataFrame:
    audit = (
        predictions.groupby(["disease", "left_out_cohort", "y_true"], observed=True)["predicted_probability"]
        .agg(
            n="size",
            mean="mean",
            median="median",
            q25=lambda x: x.quantile(0.25),
            q75=lambda x: x.quantile(0.75),
            min="min",
            max="max",
        )
        .reset_index()
    )
    audit["class"] = np.where(audit["y_true"].eq(1), "case", "control")
    return audit[["disease", "left_out_cohort", "class", "y_true", "n", "mean", "median", "q25", "q75", "min", "max"]]


def write_legend(auc_plot: pd.DataFrame, counts: pd.DataFrame, overlap: pd.DataFrame) -> None:
    med = auc_plot.groupby("disease", observed=True)["auc"].median().to_dict()
    final_k = auc_plot.groupby("disease", observed=True)["final_k"].first().to_dict()
    n_shared_2 = int(counts.loc[counts["n_diseases"].ge(2), "species_count"].sum())
    n_shared_3 = int(counts.loc[counts["n_diseases"].ge(3), "species_count"].sum())
    n_shared_4 = int(counts.loc[counts["n_diseases"].eq(4), "species_count"].sum())
    four = overlap[overlap["diseases"].eq("CRA|CRC|IBD|IBS")].iloc[0]
    legend = f"""# Feature Panel AUC and IBS-Related Overlap

**Panel A.** Locked balanced species panels were evaluated in held-out cohorts. Each row reports the held-out AUC with bootstrap 95% CI and the class-stratified predicted-probability median/IQR (control in blue-gray, case in peach); dashed lines mark AUC or probability 0.5. Final panel sizes were CRA k={int(final_k['CRA'])}, CRC k={int(final_k['CRC'])} and IBD k={int(final_k['IBD'])}.

**Panel B.** Pairwise enrichment among disease-associated feature sets. Cells show shared species count and observed/expected enrichment.

**Panel C.** Three-disease enrichment among disease-associated feature sets. Four-disease shared species are omitted from the plotted panel as requested.

Note: direct overlap between the strict final locked panels and IBS-related species is reported separately in `feature_panel_locked_panel_vs_ibs_overlap_audit.tsv`; the broad overlap claim is supported by the disease-associated feature-set analysis rather than by strict locked-panel membership alone.
"""
    (OUTDIR / "feature_panel_auc_overlap_story_legend.md").write_text(legend)


def main() -> None:
    setup_style()
    auc, predictions, final_features, assoc, overlap = load_inputs()
    sets = associated_sets(assoc)
    counts = membership_counts(sets)

    fig = plt.figure(figsize=(18.0, 5.0))
    outer = GridSpec(1, 2, figure=fig, width_ratios=[2.43, 0.94], wspace=0.085)
    auc_plot = draw_auc_forest(fig, outer[0], auc, predictions)
    pairwise, multi = draw_overlap_panel(fig, outer[1], counts, overlap, sets)
    fig.subplots_adjust(left=0.062, right=0.986, top=0.885, bottom=0.17)

    fig.savefig(OUTDIR / "feature_panel_auc_overlap_story.pdf")
    fig.savefig(OUTDIR / "feature_panel_auc_overlap_story.svg")
    fig.savefig(OUTDIR / "feature_panel_auc_overlap_story.png", dpi=300)
    fig.savefig(
        OUTDIR / "feature_panel_auc_overlap_story.tiff",
        dpi=600,
        pil_kwargs={"compression": "tiff_lzw"},
    )
    plt.close(fig)

    auc_plot.drop(columns=["plot_y"], errors="ignore").to_csv(OUTDIR / "feature_panel_auc_overlap_story_auc_source.tsv", sep="\t", index=False)
    counts.to_csv(OUTDIR / "feature_panel_auc_overlap_story_membership_counts.tsv", sep="\t", index=False)
    pairwise.to_csv(OUTDIR / "feature_panel_auc_overlap_story_pairwise_overlap.tsv", sep="\t", index=False)
    multi.to_csv(OUTDIR / "feature_panel_auc_overlap_story_multiway_overlap.tsv", sep="\t", index=False)
    locked_panel_overlap_audit(final_features, assoc).to_csv(
        OUTDIR / "feature_panel_locked_panel_vs_ibs_overlap_audit.tsv",
        sep="\t",
        index=False,
    )
    prediction_distribution_audit(predictions).to_csv(
        OUTDIR / "feature_panel_auc_overlap_story_prediction_distribution_audit.tsv",
        sep="\t",
        index=False,
    )
    write_legend(auc_plot, counts, overlap)


if __name__ == "__main__":
    main()
