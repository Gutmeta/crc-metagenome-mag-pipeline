#!/usr/bin/env python3
"""Generate deterministic, entirely synthetic inputs for every plotting example."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


SEED = 42


def write_table(frame: pd.DataFrame, path: Path) -> None:
    frame.to_csv(path, sep="\t", index=False)


def network_inputs(output_dir: Path, rng: np.random.Generator) -> None:
    count = 14
    angles = np.linspace(0, 2 * np.pi, count, endpoint=False)
    nodes = pd.DataFrame(
        {
            "node_id": [f"demo_{index:03d}" for index in range(1, count + 1)],
            "x": np.cos(angles),
            "y": np.sin(angles),
            "group": ["Guild A" if index < count / 2 else "Guild B" for index in range(count)],
            "size": rng.integers(20, 100, count),
        }
    )
    edges = []
    for index in range(count):
        for step in (1, 3):
            target = (index + step) % count
            weight = float(rng.uniform(0.25, 0.85))
            sign = "positive" if (index + step) % 3 else "negative"
            edges.append(
                {
                    "source": nodes.loc[index, "node_id"],
                    "target": nodes.loc[target, "node_id"],
                    "weight": weight,
                    "sign": sign,
                    "support": int(rng.integers(3, 9)),
                }
            )
    write_table(nodes, output_dir / "network_nodes.tsv")
    write_table(pd.DataFrame(edges).drop_duplicates(["source", "target"]), output_dir / "network_edges.tsv")


def association_inputs(output_dir: Path, rng: np.random.Generator) -> None:
    rows = ["C1alpha fraction", "C1beta fraction", "Guild balance", "Total abundance"]
    columns = ["Inflammation", "Severity", "Biomarker"]
    heatmap = []
    for row_index, row in enumerate(rows):
        for column_index, column in enumerate(columns):
            effect = float(np.clip(rng.normal((row_index - column_index) * 0.08, 0.25), -0.8, 0.8))
            heatmap.append(
                {
                    "row": row,
                    "column": column,
                    "effect": effect,
                    "q_value": float(min(0.95, np.exp(-abs(effect) * 8))),
                }
            )
    write_table(pd.DataFrame(heatmap), output_dir / "association_heatmap.tsv")

    scatter = []
    for group_index, group in enumerate(("Cohort A", "Cohort B", "Cohort C")):
        x = rng.uniform(0.0, 1.0, 24)
        y = 0.25 + (0.55 - group_index * 0.12) * x + rng.normal(0, 0.11, len(x))
        for sample_index, (x_value, y_value) in enumerate(zip(x, y), start=1):
            scatter.append(
                {
                    "sample_id": f"demo_{group_index + 1}_{sample_index:03d}",
                    "group": group,
                    "x": x_value,
                    "y": y_value,
                }
            )
    write_table(pd.DataFrame(scatter), output_dir / "association_scatter.tsv")


def performance_and_forest_inputs(output_dir: Path, rng: np.random.Generator) -> None:
    datasets = [f"Demo cohort {letter}" for letter in "ABCDEF"]
    model_a = rng.uniform(0.66, 0.88, len(datasets))
    model_b = np.clip(model_a + rng.normal(0.035, 0.045, len(datasets)), 0.55, 0.97)
    write_table(
        pd.DataFrame({"dataset": datasets, "model_a_auc": model_a, "model_b_auc": model_b}),
        output_dir / "performance_matrix.tsv",
    )

    estimates = rng.uniform(0.68, 0.91, len(datasets))
    half_width = rng.uniform(0.035, 0.08, len(datasets))
    forest = pd.DataFrame(
        {
            "comparison": [f"Comparison {index}" for index in range(1, len(datasets) + 1)],
            "estimate": estimates,
            "ci_low": np.maximum(0.5, estimates - half_width),
            "ci_high": np.minimum(0.99, estimates + half_width),
            "n_positive": rng.integers(25, 80, len(datasets)),
            "n_negative": rng.integers(25, 80, len(datasets)),
        }
    )
    write_table(forest, output_dir / "auc_forest.tsv")


def roc_inputs(output_dir: Path, rng: np.random.Generator) -> None:
    records = []
    for comparison_index, comparison in enumerate(("Model A", "Model B", "Model C")):
        labels = np.repeat([0, 1], 40)
        separation = 0.9 + comparison_index * 0.25
        scores = 1 / (1 + np.exp(-(rng.normal((labels * 2 - 1) * separation, 1.0))))
        for item_index, (label, score) in enumerate(zip(labels, scores), start=1):
            records.append(
                {
                    "sample_id": f"demo_{comparison_index + 1}_{item_index:03d}",
                    "comparison": comparison,
                    "label": int(label),
                    "score": float(score),
                }
            )
    write_table(pd.DataFrame(records), output_dir / "roc_predictions.tsv")


def structure_function_inputs(output_dir: Path, rng: np.random.Generator) -> None:
    sets = ("Set A", "Set B", "Set C", "Set D")
    memberships = []
    for item_index in range(1, 41):
        selected = [name for name in sets if rng.random() < 0.45]
        if not selected:
            selected = [sets[item_index % len(sets)]]
        for set_name in selected:
            memberships.append({"item_id": f"demo_{item_index:03d}", "set_name": set_name})
    write_table(pd.DataFrame(memberships), output_dir / "set_membership.tsv")

    sample_ids = [f"demo_{index:03d}" for index in range(1, 25)]
    groups = np.repeat(["Guild A", "Guild B", "Mixed"], 8)
    centers = {"Guild A": (-1.0, 0.2), "Guild B": (1.0, 0.0), "Mixed": (0.0, 1.0)}
    coordinates = np.array([np.array(centers[group]) + rng.normal(0, 0.35, 2) for group in groups])
    distances = np.linalg.norm(coordinates[:, None, :] - coordinates[None, :, :], axis=2)
    distance_frame = pd.DataFrame(distances, columns=sample_ids)
    distance_frame.insert(0, "sample_id", sample_ids)
    write_table(distance_frame, output_dir / "distance_matrix.tsv")
    write_table(pd.DataFrame({"sample_id": sample_ids, "group": groups}), output_dir / "ordination_metadata.tsv")

    bubble = []
    for database in ("KO", "Pathway"):
        for feature_index in range(1, 6):
            for comparison in ("Disease vs control", "Stage 2 vs stage 1", "External validation"):
                effect = float(rng.normal(0, 0.7))
                bubble.append(
                    {
                        "feature": f"{database} feature {feature_index}",
                        "database": database,
                        "comparison": comparison,
                        "effect": effect,
                        "q_value": float(min(0.99, np.exp(-abs(effect) * 5))),
                    }
                )
    write_table(pd.DataFrame(bubble), output_dir / "evidence_bubbles.tsv")


def r_plot_inputs(output_dir: Path, rng: np.random.Generator) -> None:
    tips = [f"demo_{index:03d}" for index in range(1, 13)]
    clade_a = ",".join(f"{tip}:0.2" for tip in tips[:6])
    clade_b = ",".join(f"{tip}:0.2" for tip in tips[6:])
    (output_dir / "demo_tree.nwk").write_text(f"(({clade_a}):0.3,({clade_b}):0.3);\n", encoding="utf-8")
    annotation = pd.DataFrame(
        {
            "tip_id": tips,
            "group": ["Guild A"] * 6 + ["Guild B"] * 6,
            "status": ["Cultured", "Uncultured"] * 6,
            "score": rng.uniform(0.1, 1.0, len(tips)),
        }
    )
    write_table(annotation, output_dir / "tree_annotation.tsv")

    distributions = []
    for cohort_index, cohort in enumerate(("Cohort A", "Cohort B", "Cohort C")):
        for metric_index, metric in enumerate(("C1alpha fraction", "Guild balance")):
            for group_index, group in enumerate(("Control", "Case")):
                values = rng.normal(0.35 + metric_index * 0.3 + group_index * 0.16 + cohort_index * 0.03, 0.1, 22)
                for sample_index, value in enumerate(values, start=1):
                    distributions.append(
                        {
                            "sample_id": f"demo_{cohort_index + 1}_{metric_index + 1}_{group_index + 1}_{sample_index:03d}",
                            "cohort": cohort,
                            "metric": metric,
                            "group": group,
                            "value": float(value),
                        }
                    )
    write_table(pd.DataFrame(distributions), output_dir / "group_distributions.tsv")

    volcano = pd.DataFrame(
        {
            "feature": [f"demo_feature_{index:03d}" for index in range(1, 81)],
            "effect": rng.normal(0, 0.8, 80),
            "q_value": np.clip(rng.beta(0.7, 4.0, 80), 0.0001, 0.99),
        }
    )
    volcano["label"] = np.where(volcano["q_value"] < 0.015, volcano["feature"], "")
    write_table(volcano, output_dir / "volcano.tsv")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(SEED)
    network_inputs(args.output_dir, rng)
    association_inputs(args.output_dir, rng)
    performance_and_forest_inputs(args.output_dir, rng)
    roc_inputs(args.output_dir, rng)
    structure_function_inputs(args.output_dir, rng)
    r_plot_inputs(args.output_dir, rng)
    print(f"Synthetic plotting inputs written to {args.output_dir}")


if __name__ == "__main__":
    main()
