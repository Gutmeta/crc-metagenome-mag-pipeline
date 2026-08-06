#!/usr/bin/env python3
"""Shared helpers for the portable plotting examples."""

from __future__ import annotations

import os
import random
import tempfile
from pathlib import Path
from typing import Iterable, Sequence

os.environ.setdefault(
    "MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "crc-mag-plotting-matplotlib")
)

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


DEFAULT_FORMATS = ("pdf", "svg", "png")
ALLOWED_FORMATS = set(DEFAULT_FORMATS)


def configure_style(seed: int = 42) -> None:
    """Apply deterministic, editable-vector-friendly plotting defaults."""
    random.seed(seed)
    np.random.seed(seed)
    plt.rcParams.update(
        {
            "font.family": "sans-serif",
            "font.sans-serif": ["Arial", "Liberation Sans", "DejaVu Sans"],
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.unicode_minus": False,
            "figure.facecolor": "white",
            "savefig.facecolor": "white",
        }
    )


def read_table(path: Path | str) -> pd.DataFrame:
    """Read a CSV or TSV table based on its filename suffix."""
    table_path = Path(path)
    if not table_path.is_file():
        raise FileNotFoundError(f"Input table does not exist: {table_path}")
    separator = "," if table_path.suffix.lower() == ".csv" else "\t"
    frame = pd.read_csv(table_path, sep=separator)
    if frame.empty:
        raise ValueError(f"Input table is empty: {table_path}")
    return frame


def require_columns(frame: pd.DataFrame, columns: Iterable[str], table_name: str) -> None:
    """Fail with a concise schema error if required columns are absent."""
    required = list(columns)
    missing = [column for column in required if column not in frame.columns]
    if missing:
        raise ValueError(f"{table_name} is missing required columns: {', '.join(missing)}")


def finite_numeric(frame: pd.DataFrame, columns: Iterable[str], table_name: str) -> None:
    """Coerce selected columns to numeric and reject missing or infinite values."""
    columns = list(columns)
    for column in columns:
        frame[column] = pd.to_numeric(frame[column], errors="coerce")
    if not np.isfinite(frame[columns].to_numpy(dtype=float)).all():
        raise ValueError(f"{table_name} contains missing or non-finite numeric values")


def nonempty_text(frame: pd.DataFrame, columns: Iterable[str], table_name: str) -> None:
    """Normalize selected text columns and reject missing or blank values."""
    columns = list(columns)
    for column in columns:
        if frame[column].isna().any():
            raise ValueError(f"{table_name}.{column} must not contain missing values")
        frame[column] = frame[column].astype(str)
        if frame[column].str.strip().eq("").any():
            raise ValueError(f"{table_name}.{column} must not contain empty values")


def parse_formats(value: str | Sequence[str]) -> tuple[str, ...]:
    """Parse and validate a comma-separated output format list."""
    if isinstance(value, str):
        formats = [item.strip().lower() for item in value.split(",") if item.strip()]
    else:
        formats = [str(item).strip().lower() for item in value if str(item).strip()]
    unique = tuple(dict.fromkeys(formats))
    if not unique:
        raise ValueError("At least one output format is required")
    unknown = sorted(set(unique) - ALLOWED_FORMATS)
    if unknown:
        raise ValueError(f"Unsupported output formats: {', '.join(unknown)}")
    return unique


def save_figure(
    figure: plt.Figure,
    output_prefix: Path | str,
    formats: str | Sequence[str] = DEFAULT_FORMATS,
    *,
    dpi: int = 300,
) -> list[Path]:
    """Save one figure to the requested editable vector and raster formats."""
    prefix = Path(output_prefix)
    if prefix.suffix.lower().lstrip(".") in ALLOWED_FORMATS:
        prefix = prefix.with_suffix("")
    prefix.parent.mkdir(parents=True, exist_ok=True)
    outputs: list[Path] = []
    for output_format in parse_formats(formats):
        output_path = prefix.with_suffix(f".{output_format}")
        figure.savefig(output_path, bbox_inches="tight", dpi=dpi)
        if not output_path.is_file() or output_path.stat().st_size == 0:
            raise RuntimeError(f"Figure was not written correctly: {output_path}")
        outputs.append(output_path)
    plt.close(figure)
    return outputs


def ordered_unique(values: Iterable[object]) -> list[str]:
    """Return non-empty strings in first-seen order."""
    return list(dict.fromkeys(str(value) for value in values if str(value).strip()))


def categorical_colors(values: Iterable[object]) -> dict[str, tuple[float, float, float, float]]:
    """Create a stable categorical color mapping."""
    categories = ordered_unique(values)
    color_map = plt.get_cmap("tab10")
    return {category: color_map(index % color_map.N) for index, category in enumerate(categories)}
