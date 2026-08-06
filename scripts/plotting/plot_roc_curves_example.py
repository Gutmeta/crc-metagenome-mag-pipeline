#!/usr/bin/env python3
"""Plot one or more ROC curves from out-of-fold or external predictions.

Representative source lineage: TCG-FDC editable-vector ROC figures.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import plotting_common  # Configure a writable Matplotlib cache before pyplot import.
import matplotlib.pyplot as plt
import numpy as np
from sklearn.metrics import roc_auc_score, roc_curve

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


def stratified_bootstrap_auc(
    labels: np.ndarray,
    scores: np.ndarray,
    repetitions: int,
    rng: np.random.Generator,
) -> tuple[float, float]:
    negative = np.flatnonzero(labels == 0)
    positive = np.flatnonzero(labels == 1)
    estimates = np.empty(repetitions, dtype=float)
    for repetition in range(repetitions):
        indices = np.concatenate(
            [
                rng.choice(negative, size=len(negative), replace=True),
                rng.choice(positive, size=len(positive), replace=True),
            ]
        )
        estimates[repetition] = roc_auc_score(labels[indices], scores[indices])
    return tuple(np.quantile(estimates, [0.025, 0.975]))


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, required=True)
    parser.add_argument("--output-prefix", type=Path, required=True)
    parser.add_argument("--title", default="Receiver operating characteristic")
    parser.add_argument("--formats", default="pdf,svg,png")
    parser.add_argument("--bootstrap", type=int, default=500)
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()
    if args.bootstrap < 1:
        raise ValueError("--bootstrap must be at least 1")

    configure_style(args.seed)
    predictions = read_table(args.input)
    require_columns(predictions, ["comparison", "label", "score"], "ROC predictions")
    finite_numeric(predictions, ["label", "score"], "ROC predictions")
    nonempty_text(predictions, ["comparison"], "ROC predictions")
    if not predictions["label"].isin([0, 1]).all():
        raise ValueError("label must contain only 0 and 1")
    predictions["label"] = predictions["label"].astype(int)
    if ((predictions["score"] < 0) | (predictions["score"] > 1)).any():
        raise ValueError("score must be between 0 and 1")

    colors = categorical_colors(predictions["comparison"])
    rng = np.random.default_rng(args.seed)
    figure, axis = plt.subplots(figsize=(6.2, 6.0))
    for comparison, group in predictions.groupby("comparison", sort=False):
        if set(group["label"]) != {0, 1}:
            raise ValueError(f"Comparison {comparison!r} must contain both outcome classes")
        labels = group["label"].to_numpy(dtype=int)
        scores = group["score"].to_numpy(dtype=float)
        false_positive_rate, true_positive_rate, _ = roc_curve(labels, scores)
        estimate = roc_auc_score(labels, scores)
        ci_low, ci_high = stratified_bootstrap_auc(labels, scores, args.bootstrap, rng)
        axis.plot(
            false_positive_rate,
            true_positive_rate,
            linewidth=2.2,
            color=colors[str(comparison)],
            label=f"{comparison}: {estimate:.3f} [{ci_low:.3f}, {ci_high:.3f}]",
        )

    axis.plot([0, 1], [0, 1], color="#9A9A9A", linestyle="--", linewidth=1)
    axis.set_xlim(0, 1)
    axis.set_ylim(0, 1.02)
    axis.set_aspect("equal", adjustable="box")
    axis.set_xlabel("False-positive rate")
    axis.set_ylabel("True-positive rate")
    axis.set_title(args.title, fontsize=13, weight="bold")
    axis.grid(color="#E7E7E7", linewidth=0.7)
    axis.legend(loc="lower right", frameon=False, fontsize=8.5, title="AUROC [95% CI]")
    save_figure(figure, args.output_prefix, parse_formats(args.formats))


if __name__ == "__main__":
    main()
