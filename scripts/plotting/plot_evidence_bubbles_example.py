#!/usr/bin/env python3
"""Plot functional evidence as effect-colored, FDR-sized bubbles.

Analysis context: functional evidence synthesis.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import plotting_common  # Configure a writable Matplotlib cache before pyplot import.
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import TwoSlopeNorm

from plotting_common import (
    configure_style,
    finite_numeric,
    nonempty_text,
    ordered_unique,
    parse_formats,
    read_table,
    require_columns,
    save_figure,
)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output-prefix", type=Path, required=True)
    parser.add_argument("--title", default="Functional evidence summary")
    parser.add_argument("--formats", default="pdf,svg,png")
    args = parser.parse_args()

    configure_style()
    data = read_table(args.input)
    require_columns(data, ["feature", "database", "comparison", "effect", "q_value"], "evidence bubbles")
    finite_numeric(data, ["effect", "q_value"], "evidence bubbles")
    nonempty_text(data, ["feature", "database", "comparison"], "evidence bubbles")
    if data.duplicated(["feature", "database", "comparison"]).any():
        raise ValueError("evidence bubbles contains duplicate feature/database/comparison rows")
    if ((data["q_value"] <= 0) | (data["q_value"] > 1)).any():
        raise ValueError("q_value must be greater than 0 and no greater than 1")

    comparisons = ordered_unique(data["comparison"])
    row_keys = ordered_unique(data.apply(lambda row: f"{row['database']} | {row['feature']}", axis=1))
    row_index = {key: index for index, key in enumerate(row_keys)}
    comparison_index = {key: index for index, key in enumerate(comparisons)}
    y = np.array([row_index[f"{row.database} | {row.feature}"] for row in data.itertuples(index=False)])
    x = np.array([comparison_index[str(value)] for value in data["comparison"]])
    significance = -np.log10(data["q_value"].to_numpy(dtype=float))
    bubble_sizes = 25 + 150 * significance / max(1.0, float(significance.max()))
    effect_limit = max(0.1, float(data["effect"].abs().max()))

    figure, axis = plt.subplots(figsize=(max(7.0, 1.5 * len(comparisons) + 3.5), max(5.0, 0.38 * len(row_keys) + 1.8)))
    scatter = axis.scatter(
        x,
        y,
        s=bubble_sizes,
        c=data["effect"],
        cmap="RdBu_r",
        norm=TwoSlopeNorm(vmin=-effect_limit, vcenter=0, vmax=effect_limit),
        edgecolor="#404040",
        linewidth=0.4,
        alpha=0.88,
    )
    axis.set_xticks(range(len(comparisons)), comparisons, rotation=30, ha="right")
    axis.set_yticks(range(len(row_keys)), row_keys)
    axis.set_ylim(len(row_keys) - 0.5, -0.5)
    axis.set_title(args.title, fontsize=13, weight="bold")
    axis.grid(color="#E8E8E8", linewidth=0.65)
    axis.tick_params(length=0)
    colorbar = figure.colorbar(scatter, ax=axis, pad=0.02, shrink=0.78)
    colorbar.set_label("Effect size")
    legend_values = [1, 2, 3]
    handles = [
        axis.scatter([], [], s=25 + 150 * value / max(1.0, float(significance.max())), color="#888888", alpha=0.75)
        for value in legend_values
    ]
    axis.legend(handles, [f"{value}" for value in legend_values], title="-log10(q)", frameon=False, loc="upper left", bbox_to_anchor=(1.18, 1.0))
    save_figure(figure, args.output_prefix, parse_formats(args.formats))


if __name__ == "__main__":
    main()
