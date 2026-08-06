#!/usr/bin/env python3
"""Plot grouped associations with linear fits and confidence bands.

Representative source lineage: pipeline_IBD within-cohort association figures.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import plotting_common  # Configure a writable Matplotlib cache before pyplot import.
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import linregress

from plotting_common import (
    categorical_colors,
    configure_style,
    finite_numeric,
    parse_formats,
    read_table,
    require_columns,
    save_figure,
)


def regression_band(x: np.ndarray, y: np.ndarray, grid: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    fit = linregress(x, y)
    predicted = fit.intercept + fit.slope * grid
    residuals = y - (fit.intercept + fit.slope * x)
    residual_se = np.sqrt(np.sum(residuals**2) / (len(x) - 2))
    denominator = np.sum((x - np.mean(x)) ** 2)
    if denominator <= 0:
        raise ValueError("Each group must contain variation in x")
    confidence = 1.96 * residual_se * np.sqrt(1 / len(x) + (grid - np.mean(x)) ** 2 / denominator)
    return predicted, confidence


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output-prefix", type=Path, required=True)
    parser.add_argument("--title", default="Within-cohort association")
    parser.add_argument("--x-label", default="Predictor")
    parser.add_argument("--y-label", default="Outcome")
    parser.add_argument("--formats", default="pdf,svg,png")
    args = parser.parse_args()

    configure_style()
    data = read_table(args.input)
    require_columns(data, ["x", "y", "group"], "association scatter")
    finite_numeric(data, ["x", "y"], "association scatter")
    if data["group"].astype(str).str.strip().eq("").any():
        raise ValueError("group must not contain empty values")

    colors = categorical_colors(data["group"])
    figure, axis = plt.subplots(figsize=(7.2, 5.5))
    for group, group_frame in data.groupby("group", sort=False):
        if len(group_frame) < 4:
            raise ValueError(f"Group {group!r} requires at least four observations")
        x = group_frame["x"].to_numpy(dtype=float)
        y = group_frame["y"].to_numpy(dtype=float)
        grid = np.linspace(float(x.min()), float(x.max()), 100)
        predicted, confidence = regression_band(x, y, grid)
        fit = linregress(x, y)
        color = colors[str(group)]
        axis.scatter(x, y, s=30, alpha=0.72, color=color, edgecolor="white", linewidth=0.4)
        axis.plot(grid, predicted, color=color, lw=2.0, label=f"{group}: r={fit.rvalue:.2f}, p={fit.pvalue:.3g}")
        axis.fill_between(grid, predicted - confidence, predicted + confidence, color=color, alpha=0.14, linewidth=0)

    axis.axhline(0, color="#B8B8B8", lw=0.7, zorder=0)
    axis.set_xlabel(args.x_label)
    axis.set_ylabel(args.y_label)
    axis.set_title(args.title, fontsize=13, weight="bold")
    axis.grid(color="#E6E6E6", linewidth=0.6, alpha=0.8)
    axis.legend(frameon=False, fontsize=9)
    save_figure(figure, args.output_prefix, parse_formats(args.formats))


if __name__ == "__main__":
    main()
