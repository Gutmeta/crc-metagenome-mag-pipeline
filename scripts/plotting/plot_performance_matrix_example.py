#!/usr/bin/env python3
"""Plot a compact annotated matrix comparing two model AUROCs.

Representative source lineage: pipeline_Pan performance-comparison figures.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import plotting_common  # Configure a writable Matplotlib cache before pyplot import.
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import Normalize, TwoSlopeNorm

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
    parser.add_argument("--title", default="Model performance comparison")
    parser.add_argument("--model-a-label", default="Model A AUROC")
    parser.add_argument("--model-b-label", default="Model B AUROC")
    parser.add_argument("--formats", default="pdf,svg,png")
    args = parser.parse_args()

    configure_style()
    data = read_table(args.input)
    require_columns(data, ["dataset", "model_a_auc", "model_b_auc"], "performance matrix")
    finite_numeric(data, ["model_a_auc", "model_b_auc"], "performance matrix")
    nonempty_text(data, ["dataset"], "performance matrix")
    if data["dataset"].duplicated().any():
        raise ValueError("performance matrix contains duplicate dataset values")
    values = data[["model_a_auc", "model_b_auc"]].to_numpy(dtype=float)
    if ((values < 0) | (values > 1)).any():
        raise ValueError("AUROC values must be between 0 and 1")
    data["delta"] = data["model_b_auc"] - data["model_a_auc"]

    columns = [
        ("model_a_auc", args.model_a_label, plt.get_cmap("Blues"), Normalize(vmin=0.5, vmax=1.0)),
        ("model_b_auc", args.model_b_label, plt.get_cmap("Purples"), Normalize(vmin=0.5, vmax=1.0)),
    ]
    delta_limit = max(0.02, float(data["delta"].abs().max()))
    columns.append(("delta", "Delta (B - A)", plt.get_cmap("RdBu_r"), TwoSlopeNorm(vmin=-delta_limit, vcenter=0, vmax=delta_limit)))

    figure, axes = plt.subplots(1, 3, figsize=(8.2, max(3.6, 0.58 * len(data) + 1.6)), sharey=True)
    for axis_index, (column, label, color_map, norm) in enumerate(columns):
        axis = axes[axis_index]
        column_values = data[column].to_numpy(dtype=float)[:, None]
        axis.imshow(column_values, cmap=color_map, norm=norm, aspect="auto")
        axis.set_xticks([0], [label], rotation=25, ha="right")
        axis.set_yticks(range(len(data)))
        if axis_index == 0:
            axis.set_yticklabels(data["dataset"])
        axis.tick_params(length=0)
        for row_index, value in enumerate(data[column]):
            label_text = f"{value:+.3f}" if column == "delta" else f"{value:.3f}"
            axis.text(0, row_index, label_text, ha="center", va="center", fontsize=9, color="#181818")
        for spine in axis.spines.values():
            spine.set_visible(False)
    figure.suptitle(args.title, fontsize=13, weight="bold", y=1.01)
    figure.subplots_adjust(wspace=0.08)
    save_figure(figure, args.output_prefix, parse_formats(args.formats))


if __name__ == "__main__":
    main()
