#!/usr/bin/env python3
"""Plot an UpSet-style intersection summary from long-form set membership.

Representative source lineage: TCG-FDC structure/function overlap figures.
"""

from __future__ import annotations

import argparse
from collections import Counter
from pathlib import Path

import plotting_common  # Configure a writable Matplotlib cache before pyplot import.
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.gridspec import GridSpec

from plotting_common import (
    configure_style,
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
    parser.add_argument("--title", default="Set intersections")
    parser.add_argument("--formats", default="pdf,svg,png")
    parser.add_argument("--max-intersections", type=int, default=15)
    args = parser.parse_args()
    if args.max_intersections < 1:
        raise ValueError("--max-intersections must be positive")

    configure_style()
    membership = read_table(args.input)
    require_columns(membership, ["item_id", "set_name"], "set membership")
    membership = membership[["item_id", "set_name"]].copy()
    nonempty_text(membership, ["item_id", "set_name"], "set membership")
    membership = membership.drop_duplicates()
    set_order = ordered_unique(membership["set_name"])
    if len(set_order) < 2:
        raise ValueError("At least two sets are required")

    item_sets = membership.groupby("item_id", sort=False)["set_name"].apply(lambda values: frozenset(values))
    intersections = Counter(item_sets)
    ranked = sorted(intersections.items(), key=lambda item: (-item[1], -len(item[0]), sorted(item[0])))
    ranked = ranked[: args.max_intersections]
    combinations = [combination for combination, _ in ranked]
    counts = np.array([count for _, count in ranked], dtype=int)

    figure = plt.figure(figsize=(max(7.0, 0.58 * len(ranked) + 3.0), 6.0))
    grid = GridSpec(2, 1, height_ratios=[2.0, 1.45], hspace=0.05, figure=figure)
    bar_axis = figure.add_subplot(grid[0])
    matrix_axis = figure.add_subplot(grid[1], sharex=bar_axis)
    x_positions = np.arange(len(ranked))
    bar_axis.bar(x_positions, counts, color="#4C78A8", width=0.72)
    for x_position, count in zip(x_positions, counts):
        bar_axis.text(x_position, count, str(count), ha="center", va="bottom", fontsize=8)
    bar_axis.set_ylabel("Intersection size")
    bar_axis.set_title(args.title, fontsize=13, weight="bold")
    bar_axis.grid(axis="y", color="#E5E5E5", linewidth=0.7)
    bar_axis.tick_params(axis="x", labelbottom=False, length=0)

    y_positions = np.arange(len(set_order))
    for y_position in y_positions:
        matrix_axis.axhline(y_position, color="#ECECEC", linewidth=0.7, zorder=0)
    for x_position, combination in zip(x_positions, combinations):
        active = [set_order.index(set_name) for set_name in set_order if set_name in combination]
        matrix_axis.scatter(
            np.full(len(set_order), x_position),
            y_positions,
            s=32,
            color="#D9D9D9",
            zorder=1,
        )
        matrix_axis.scatter(
            np.full(len(active), x_position),
            active,
            s=44,
            color="#202020",
            zorder=3,
        )
        if len(active) > 1:
            matrix_axis.plot([x_position, x_position], [min(active), max(active)], color="#202020", lw=1.5, zorder=2)
    matrix_axis.set_yticks(y_positions, set_order)
    matrix_axis.set_xticks(x_positions, [str(index + 1) for index in range(len(ranked))])
    matrix_axis.set_xlabel("Intersection rank")
    matrix_axis.set_ylim(len(set_order) - 0.5, -0.5)
    matrix_axis.tick_params(length=0)
    save_figure(figure, args.output_prefix, parse_formats(args.formats))


if __name__ == "__main__":
    main()
