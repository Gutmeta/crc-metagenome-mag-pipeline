#!/usr/bin/env python3
"""Build auditable cohort-level ML evidence scores for phylogeny tips."""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from plotting_common import finite_numeric, nonempty_text, read_table, require_columns


METHOD_NAME = "minimum-rank cohort evidence"
SCORE_FORMULA = "(rank_min - 1) / (N - 1); zero importance forced to zero"
COHORT_WEIGHTING = "equal"


def build_evidence(
    importance: pd.DataFrame,
    tips: pd.DataFrame,
    *,
    top_quartile_threshold: float = 0.75,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Validate inputs and return tip summaries plus the complete audit table."""
    require_columns(
        importance,
        ("feature_id", "cohort", "mean_importance"),
        "importance table",
    )
    require_columns(tips, ("tip_id",), "tree-node table")
    importance = importance.loc[:, ["feature_id", "cohort", "mean_importance"]].copy()
    tips = tips.loc[:, ["tip_id"]].copy()
    nonempty_text(importance, ("feature_id", "cohort"), "importance table")
    nonempty_text(tips, ("tip_id",), "tree-node table")
    finite_numeric(importance, ("mean_importance",), "importance table")

    if not np.isfinite(top_quartile_threshold) or not 0 <= top_quartile_threshold <= 1:
        raise ValueError("top-quartile threshold must lie in [0, 1]")
    if importance.duplicated(["feature_id", "cohort"]).any():
        raise ValueError("importance table contains duplicate feature_id/cohort keys")
    if tips["tip_id"].duplicated().any():
        raise ValueError("tree-node table contains duplicate tip_id values")
    if (importance["mean_importance"] < 0).any():
        raise ValueError("mean_importance values must be non-negative")

    features = list(dict.fromkeys(importance["feature_id"]))
    cohorts = list(dict.fromkeys(importance["cohort"]))
    if len(features) < 2:
        raise ValueError("importance table must contain at least two candidate features")
    if not cohorts:
        raise ValueError("importance table must contain at least one cohort")
    expected_keys = pd.MultiIndex.from_product(
        [features, cohorts], names=["feature_id", "cohort"]
    )
    observed_keys = pd.MultiIndex.from_frame(importance[["feature_id", "cohort"]])
    missing_keys = expected_keys.difference(observed_keys)
    if len(missing_keys):
        preview = ", ".join(f"{feature}/{cohort}" for feature, cohort in missing_keys[:5])
        raise ValueError(
            "importance table must contain a complete feature-by-cohort matrix; "
            f"missing {len(missing_keys)} key(s), including {preview}"
        )

    unknown_tips = sorted(set(tips["tip_id"]) - set(features))
    if unknown_tips:
        preview = ", ".join(unknown_tips[:5])
        raise ValueError(
            "tree-node table contains tip_id values absent from the candidate universe: "
            f"{preview}"
        )

    audit = importance.copy()
    audit["rank_min"] = audit.groupby("cohort", sort=False)["mean_importance"].rank(
        method="min", ascending=True
    ).astype(int)
    universe_n = len(features)
    audit["cohort_evidence_score"] = (audit["rank_min"] - 1.0) / (universe_n - 1.0)
    audit.loc[audit["mean_importance"].eq(0), "cohort_evidence_score"] = 0.0
    audit["top_quartile_supported"] = audit["cohort_evidence_score"].ge(
        top_quartile_threshold
    )
    audit["nonzero_supported"] = audit["mean_importance"].gt(0)
    audit["included_tree_tip"] = audit["feature_id"].isin(tips["tip_id"])
    audit["normalization_universe_n"] = universe_n
    audit["normalization_method"] = METHOD_NAME
    audit["score_formula"] = SCORE_FORMULA
    audit["cohort_weighting"] = COHORT_WEIGHTING
    audit["top_quartile_threshold"] = top_quartile_threshold
    audit = audit.sort_values(
        ["cohort", "rank_min", "feature_id"], kind="stable"
    ).reset_index(drop=True)

    tip_audit = audit[audit["included_tree_tip"]]
    summary = (
        tip_audit.groupby("feature_id", sort=False)
        .agg(
            cohort_evidence_score=("cohort_evidence_score", "mean"),
            top_quartile_support_count=("top_quartile_supported", "sum"),
            nonzero_support_count=("nonzero_supported", "sum"),
        )
        .reindex(tips["tip_id"])
        .reset_index(names="tip_id")
    )
    summary["cohort_count"] = len(cohorts)
    summary["normalization_universe_n"] = universe_n
    summary["normalization_method"] = METHOD_NAME
    summary["score_formula"] = SCORE_FORMULA
    summary["cohort_weighting"] = COHORT_WEIGHTING
    summary["top_quartile_threshold"] = top_quartile_threshold
    summary["top_quartile_support_count"] = summary[
        "top_quartile_support_count"
    ].astype(int)
    summary["nonzero_support_count"] = summary["nonzero_support_count"].astype(int)
    return summary, audit


def write_table(frame: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(path, sep="\t", index=False)
    if not path.is_file() or path.stat().st_size == 0:
        raise RuntimeError(f"Output table was not written correctly: {path}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--importance-input", type=Path, required=True)
    parser.add_argument("--tips-input", type=Path, required=True)
    parser.add_argument("--summary-output", type=Path, required=True)
    parser.add_argument("--audit-output", type=Path, required=True)
    parser.add_argument("--top-quartile-threshold", type=float, default=0.75)
    args = parser.parse_args()

    summary, audit = build_evidence(
        read_table(args.importance_input),
        read_table(args.tips_input),
        top_quartile_threshold=args.top_quartile_threshold,
    )
    write_table(summary, args.summary_output)
    write_table(audit, args.audit_output)
    print(f"Cohort-evidence summary written to {args.summary_output}")
    print(f"Complete audit table written to {args.audit_output}")


if __name__ == "__main__":
    main()
