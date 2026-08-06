#!/usr/bin/env python3
"""Plot an effect-size heatmap with FDR annotations.

Representative source lineage: pipeline_IBD clinical-association figures.
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
    ordered_unique,
    parse_formats,
    read_table,
    require_columns,
    save_figure,
)


def significance_stars(q_value: float) -> str:
    if q_value < 0.001:
        return "***"
    if q_value < 0.01:
        return "**"
    if q_value < 0.05:
        return "*"
    return ""


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output-prefix", type=Path, required=True)
    parser.add_argument("--title", default="Clinical association heatmap")
    parser.add_argument("--formats", default="pdf,svg,png")
    args = parser.parse_args()

    configure_style()
    data = read_table(args.input)
    require_columns(data, ["row", "column", "effect", "q_value"], "association heatmap")
    finite_numeric(data, ["effect", "q_value"], "association heatmap")
    if data.duplicated(["row", "column"]).any():
        raise ValueError("association heatmap contains duplicate row/column pairs")
    if ((data["q_value"] < 0) | (data["q_value"] > 1)).any():
        raise ValueError("q_value must be between 0 and 1")

    row_order = ordered_unique(data["row"])
    column_order = ordered_unique(data["column"])
    effects = data.pivot(index="row", columns="column", values="effect").reindex(index=row_order, columns=column_order)
    q_values = data.pivot(index="row", columns="column", values="q_value").reindex(index=row_order, columns=column_order)
    limit = max(0.1, float(np.nanmax(np.abs(effects.to_numpy(dtype=float)))))

    figure_width = max(5.2, 1.25 * len(column_order) + 2.6)
    figure_height = max(3.8, 0.72 * len(row_order) + 1.8)
    figure, axis = plt.subplots(figsize=(figure_width, figure_height))
    color_map = plt.get_cmap("RdBu_r").copy()
    color_map.set_bad("#EEEEEE")
    image = axis.imshow(
        np.ma.masked_invalid(effects.to_numpy(dtype=float)),
        cmap=color_map,
        norm=TwoSlopeNorm(vmin=-limit, vcenter=0, vmax=limit),
        aspect="auto",
    )
    axis.set_xticks(range(len(column_order)), column_order, rotation=35, ha="right")
    axis.set_yticks(range(len(row_order)), row_order)
    axis.set_title(args.title, fontsize=13, weight="bold", pad=12)
    for row_index in range(len(row_order)):
        for column_index in range(len(column_order)):
            effect = effects.iloc[row_index, column_index]
            q_value = q_values.iloc[row_index, column_index]
            if not np.isfinite(effect):
                continue
            stars = significance_stars(float(q_value)) if np.isfinite(q_value) else ""
            text_color = "white" if abs(float(effect)) > limit * 0.55 else "#202020"
            axis.text(
                column_index,
                row_index,
                f"{effect:+.2f}{stars}",
                ha="center",
                va="center",
                fontsize=9,
                color=text_color,
            )
    colorbar = figure.colorbar(image, ax=axis, shrink=0.82, pad=0.03)
    colorbar.set_label("Effect size")
    axis.tick_params(length=0)
    for spine in axis.spines.values():
        spine.set_visible(False)
    save_figure(figure, args.output_prefix, parse_formats(args.formats))


if __name__ == "__main__":
    main()
