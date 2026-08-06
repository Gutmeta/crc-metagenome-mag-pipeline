#!/usr/bin/env python3
"""Plot a portable consensus co-abundance network from node and edge tables.

Analysis context: microbial co-abundance network visualization.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import plotting_common  # Configure a writable Matplotlib cache before pyplot import.
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import LineCollection
from matplotlib.lines import Line2D

from plotting_common import (
    categorical_colors,
    configure_style,
    finite_numeric,
    nonempty_text,
    parse_formats,
    read_table,
    require_columns,
    save_figure,
)


def rescale(values: np.ndarray, low: float, high: float) -> np.ndarray:
    values = np.asarray(values, dtype=float)
    minimum, maximum = float(values.min()), float(values.max())
    if np.isclose(minimum, maximum):
        return np.full(values.shape, (low + high) / 2.0)
    return low + (values - minimum) * (high - low) / (maximum - minimum)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--nodes", type=Path, required=True)
    parser.add_argument("--edges", type=Path, required=True)
    parser.add_argument("--output-prefix", type=Path, required=True)
    parser.add_argument("--title", default="Consensus co-abundance network")
    parser.add_argument("--formats", default="pdf,svg,png")
    parser.add_argument("--label-nodes", action="store_true")
    args = parser.parse_args()

    configure_style()
    nodes = read_table(args.nodes)
    edges = read_table(args.edges)
    require_columns(nodes, ["node_id", "x", "y", "group", "size"], "nodes")
    require_columns(edges, ["source", "target", "weight", "sign", "support"], "edges")
    finite_numeric(nodes, ["x", "y", "size"], "nodes")
    finite_numeric(edges, ["weight", "support"], "edges")
    nonempty_text(nodes, ["node_id", "group"], "nodes")
    nonempty_text(edges, ["source", "target", "sign"], "edges")

    if nodes["node_id"].duplicated().any():
        raise ValueError("nodes contains duplicate node_id values")
    if (nodes["size"] <= 0).any():
        raise ValueError("nodes.size must be positive")
    if ((edges["weight"] <= 0) | (edges["support"] <= 0)).any():
        raise ValueError("edges.weight and edges.support must be positive")
    allowed_signs = {"positive", "negative"}
    unknown_signs = sorted(set(edges["sign"].astype(str)) - allowed_signs)
    if unknown_signs:
        raise ValueError(f"edges.sign contains unsupported values: {', '.join(unknown_signs)}")

    node_ids = set(nodes["node_id"])
    edge_ids = set(edges["source"].astype(str)) | set(edges["target"].astype(str))
    unknown_nodes = sorted(edge_ids - node_ids)
    if unknown_nodes:
        raise ValueError(f"edges reference unknown nodes: {', '.join(unknown_nodes[:8])}")

    coordinates = nodes.set_index("node_id")[["x", "y"]]
    segments = [
        [coordinates.loc[str(row.source)].to_numpy(), coordinates.loc[str(row.target)].to_numpy()]
        for row in edges.itertuples(index=False)
    ]
    weights = np.abs(edges["weight"].to_numpy(dtype=float))
    supports = edges["support"].to_numpy(dtype=float)
    widths = rescale(weights, 0.4, 2.6)
    alphas = rescale(supports, 0.20, 0.80)
    edge_palette = {"positive": "#A85B50", "negative": "#3B78A8"}
    edge_colors = [
        (*plt.matplotlib.colors.to_rgb(edge_palette[str(sign)]), float(alpha))
        for sign, alpha in zip(edges["sign"], alphas)
    ]

    figure, axis = plt.subplots(figsize=(8.0, 7.2))
    axis.add_collection(LineCollection(segments, colors=edge_colors, linewidths=widths, zorder=1))

    group_colors = categorical_colors(nodes["group"])
    node_sizes = rescale(np.sqrt(nodes["size"].to_numpy(dtype=float)), 55, 330)
    for group, group_frame in nodes.groupby("group", sort=False):
        positions = group_frame.index.to_numpy()
        axis.scatter(
            group_frame["x"],
            group_frame["y"],
            s=node_sizes[positions],
            color=group_colors[str(group)],
            edgecolor="white",
            linewidth=0.8,
            label=str(group),
            zorder=3,
        )
        if args.label_nodes:
            label_column = "label" if "label" in nodes.columns else "node_id"
            for row in group_frame.itertuples():
                axis.text(row.x, row.y, str(getattr(row, label_column)), fontsize=7, ha="center", va="bottom")

    axis.autoscale()
    axis.margins(0.12)
    axis.set_aspect("equal", adjustable="datalim")
    axis.set_axis_off()
    axis.set_title(args.title, fontsize=14, weight="bold", pad=16)
    handles = [
        Line2D([0], [0], marker="o", linestyle="", color=color, label=group, markersize=8)
        for group, color in group_colors.items()
    ]
    handles.extend(
        [
            Line2D([0], [0], color=edge_palette["positive"], lw=2, label="Positive edge"),
            Line2D([0], [0], color=edge_palette["negative"], lw=2, label="Negative edge"),
        ]
    )
    axis.legend(handles=handles, loc="upper left", bbox_to_anchor=(1.01, 1.0), frameon=False)
    save_figure(figure, args.output_prefix, parse_formats(args.formats))


if __name__ == "__main__":
    main()
