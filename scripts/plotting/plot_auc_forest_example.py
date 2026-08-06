#!/usr/bin/env python3
"""Plot precomputed AUROC estimates and confidence intervals as a forest plot.

Analysis context: cohort-level model-discrimination summaries.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import plotting_common  # Configure a writable Matplotlib cache before pyplot import.
import matplotlib.pyplot as plt
import numpy as np

from plotting_common import (
    configure_style,
    finite_numeric,
    nonempty_text,
    parse_formats,
    read_table,
    require_columns,
    save_figure,
)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output-prefix", type=Path, required=True)
    parser.add_argument("--title", default="AUROC summary")
    parser.add_argument("--formats", default="pdf,svg,png")
    args = parser.parse_args()

    configure_style()
    data = read_table(args.input)
    required = ["comparison", "estimate", "ci_low", "ci_high", "n_positive", "n_negative"]
    require_columns(data, required, "AUC forest")
    finite_numeric(data, ["estimate", "ci_low", "ci_high", "n_positive", "n_negative"], "AUC forest")
    nonempty_text(data, ["comparison"], "AUC forest")
    if data["comparison"].duplicated().any():
        raise ValueError("AUC forest contains duplicate comparison values")
    invalid = (
        (data["ci_low"] < 0)
        | (data["ci_high"] > 1)
        | (data["ci_low"] > data["estimate"])
        | (data["estimate"] > data["ci_high"])
    )
    if invalid.any():
        raise ValueError("Each row must satisfy 0 <= ci_low <= estimate <= ci_high <= 1")
    if ((data["n_positive"] <= 0) | (data["n_negative"] <= 0)).any():
        raise ValueError("n_positive and n_negative must be positive")
    sample_counts = data[["n_positive", "n_negative"]].to_numpy(dtype=float)
    if not np.allclose(sample_counts, np.round(sample_counts)):
        raise ValueError("n_positive and n_negative must be integers")

    plot_data = data.iloc[::-1].reset_index(drop=True)
    y_positions = np.arange(len(plot_data))
    figure, axis = plt.subplots(figsize=(8.4, max(4.0, 0.62 * len(plot_data) + 1.8)))
    for y_position in y_positions[::2]:
        axis.axhspan(y_position - 0.48, y_position + 0.48, color="#F5F3F7", zorder=0)
    lower_error = plot_data["estimate"] - plot_data["ci_low"]
    upper_error = plot_data["ci_high"] - plot_data["estimate"]
    axis.errorbar(
        plot_data["estimate"],
        y_positions,
        xerr=np.vstack([lower_error, upper_error]),
        fmt="o",
        color="#6F4A8E",
        ecolor="#6F4A8E",
        markersize=6,
        capsize=3,
        linewidth=1.6,
        zorder=3,
    )
    axis.axvline(0.5, color="#9A9A9A", linestyle="--", linewidth=1)
    axis.set_yticks(y_positions, plot_data["comparison"])
    axis.set_xlim(0.45, 1.08)
    axis.set_xlabel("AUROC (95% confidence interval)")
    axis.set_title(args.title, fontsize=13, weight="bold", pad=12)
    axis.grid(axis="x", color="#E1E1E1", linewidth=0.7)
    for y_position, row in enumerate(plot_data.itertuples(index=False)):
        axis.text(
            1.075,
            y_position,
            f"{row.estimate:.3f} [{row.ci_low:.3f}, {row.ci_high:.3f}]\n"
            f"n={int(row.n_positive)}+/{int(row.n_negative)}-",
            ha="right",
            va="center",
            fontsize=8,
        )
    save_figure(figure, args.output_prefix, parse_formats(args.formats))


if __name__ == "__main__":
    main()
