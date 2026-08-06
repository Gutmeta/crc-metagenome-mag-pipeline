#!/usr/bin/env python3
"""Compute and plot classical PCoA from a square distance matrix.

Representative source lineage: TCG-FDC and pipeline_IBD ordination figures.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import plotting_common  # Configure a writable Matplotlib cache before pyplot import.
import matplotlib.pyplot as plt
import numpy as np

from plotting_common import (
    categorical_colors,
    configure_style,
    finite_numeric,
    parse_formats,
    read_table,
    require_columns,
    save_figure,
)


def classical_pcoa(distance: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    sample_count = distance.shape[0]
    centering = np.eye(sample_count) - np.ones((sample_count, sample_count)) / sample_count
    gram = -0.5 * centering @ (distance**2) @ centering
    eigenvalues, eigenvectors = np.linalg.eigh(gram)
    order = np.argsort(eigenvalues)[::-1]
    eigenvalues = eigenvalues[order]
    eigenvectors = eigenvectors[:, order]
    positive = eigenvalues > np.finfo(float).eps * max(1.0, float(eigenvalues[0]))
    if positive.sum() < 2:
        raise ValueError("Distance matrix must yield at least two positive PCoA axes")
    positive_values = eigenvalues[positive]
    coordinates = eigenvectors[:, positive] * np.sqrt(positive_values)
    explained = positive_values / positive_values.sum()
    return coordinates[:, :2], explained[:2]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--distance-matrix", type=Path, required=True)
    parser.add_argument("--metadata", type=Path, required=True)
    parser.add_argument("--output-prefix", type=Path, required=True)
    parser.add_argument("--title", default="Principal coordinates analysis")
    parser.add_argument("--formats", default="pdf,svg,png")
    args = parser.parse_args()

    configure_style()
    distance_frame = read_table(args.distance_matrix)
    if distance_frame.shape[1] < 3:
        raise ValueError("Distance matrix must include sample_id and at least two sample columns")
    identifier_column = distance_frame.columns[0]
    sample_ids = distance_frame[identifier_column].astype(str).tolist()
    if len(sample_ids) != len(set(sample_ids)):
        raise ValueError("Distance matrix contains duplicate row identifiers")
    matrix_columns = [str(column) for column in distance_frame.columns[1:]]
    if set(matrix_columns) != set(sample_ids) or len(matrix_columns) != len(sample_ids):
        raise ValueError("Distance matrix must be square with matching row and column identifiers")
    distance_frame.columns = [identifier_column, *matrix_columns]
    distance_frame = distance_frame.set_index(identifier_column).loc[sample_ids, sample_ids]
    finite_numeric(distance_frame, distance_frame.columns.tolist(), "distance matrix")
    distance = distance_frame.to_numpy(dtype=float)
    if (distance < 0).any():
        raise ValueError("Distance matrix cannot contain negative values")
    if not np.allclose(distance, distance.T, atol=1e-8):
        raise ValueError("Distance matrix must be symmetric")
    if not np.allclose(np.diag(distance), 0, atol=1e-8):
        raise ValueError("Distance matrix diagonal must be zero")

    metadata = read_table(args.metadata)
    require_columns(metadata, ["sample_id", "group"], "ordination metadata")
    metadata["sample_id"] = metadata["sample_id"].astype(str)
    if metadata["sample_id"].duplicated().any():
        raise ValueError("ordination metadata contains duplicate sample_id values")
    metadata = metadata.set_index("sample_id")
    missing = sorted(set(sample_ids) - set(metadata.index))
    if missing:
        raise ValueError(f"ordination metadata is missing samples: {', '.join(missing[:8])}")
    metadata = metadata.loc[sample_ids]

    coordinates, explained = classical_pcoa(distance)
    colors = categorical_colors(metadata["group"])
    figure, axis = plt.subplots(figsize=(6.4, 5.6))
    for group, indices in metadata.groupby("group", sort=False).groups.items():
        positions = np.array([sample_ids.index(sample_id) for sample_id in indices])
        axis.scatter(
            coordinates[positions, 0],
            coordinates[positions, 1],
            s=48,
            alpha=0.82,
            color=colors[str(group)],
            edgecolor="white",
            linewidth=0.5,
            label=str(group),
        )
    axis.axhline(0, color="#D8D8D8", linewidth=0.7)
    axis.axvline(0, color="#D8D8D8", linewidth=0.7)
    axis.set_xlabel(f"PCoA 1 ({explained[0] * 100:.1f}%)")
    axis.set_ylabel(f"PCoA 2 ({explained[1] * 100:.1f}%)")
    axis.set_title(args.title, fontsize=13, weight="bold")
    axis.legend(frameon=False)
    save_figure(figure, args.output_prefix, parse_formats(args.formats))


if __name__ == "__main__":
    main()
