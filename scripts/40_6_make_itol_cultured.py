#!/usr/bin/env python3
"""Call species-level culture representation from skani hits and write an iTOL track."""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import pandas as pd


DEFAULT_CULTURED_COLOR = "#8FBBD8"
DEFAULT_UNCULTURED_COLOR = "#C9E29A"
HEX_COLOR = re.compile(r"^#[0-9A-Fa-f]{6}$")
ACCESSION_PATTERN = re.compile(r"(GC[AF]_\d+\.\d+)")


def read_table(path: Path) -> pd.DataFrame:
    if not path.is_file():
        raise FileNotFoundError(f"Input table does not exist: {path}")
    separator = "," if path.suffix.lower() == ".csv" else "\t"
    frame = pd.read_csv(path, sep=separator, dtype=str, keep_default_na=False, low_memory=False)
    if frame.empty:
        raise ValueError(f"Input table is empty: {path}")
    return frame


def require_columns(frame: pd.DataFrame, columns: set[str], table_name: str) -> None:
    missing = sorted(columns.difference(frame.columns))
    if missing:
        raise ValueError(f"{table_name} is missing columns: {', '.join(missing)}")


def normalize_path(value: object) -> str:
    return str(Path(str(value)).expanduser().resolve())


def normalize_accession(value: object) -> str:
    return re.sub(r"^(RS_|GB_)", "", str(value).strip())


def accession_from_path(value: object) -> str:
    match = ACCESSION_PATTERN.search(Path(str(value)).name)
    return match.group(1) if match else ""


def species_from_taxonomy(value: object) -> str:
    species = [token.strip() for token in str(value).split(";") if token.strip().startswith("s__")]
    return species[-1] if species else ""


def validate_color(value: str, option: str) -> str:
    if not HEX_COLOR.fullmatch(value):
        raise ValueError(f"{option} must use #RRGGBB notation")
    return value.upper()


def load_query_map(path: Path) -> pd.DataFrame:
    frame = read_table(path)
    require_columns(frame, {"tip_id", "query_fasta"}, "query map")
    frame = frame.loc[:, ["tip_id", "query_fasta"]].copy()
    for column in ["tip_id", "query_fasta"]:
        frame[column] = frame[column].str.strip()
        if frame[column].eq("").any():
            raise ValueError(f"query map {column} must not contain empty values")
    if frame["tip_id"].duplicated().any():
        raise ValueError("query map contains duplicate tip_id values")
    frame["query_path_norm"] = frame["query_fasta"].map(normalize_path)
    return frame


def load_skani_results(path: Path, known_query_paths: set[str]) -> pd.DataFrame:
    frame = read_table(path)
    expected = {"Ref_file", "Query_file", "ANI", "Align_fraction_ref", "Align_fraction_query"}
    require_columns(frame, expected, "skani results")
    frame = frame.loc[:, list(expected)].rename(
        columns={
            "Ref_file": "reference_fasta",
            "Query_file": "query_fasta",
            "ANI": "ani",
            "Align_fraction_ref": "af_reference_percent",
            "Align_fraction_query": "af_query_percent",
        }
    )
    for column in ["ani", "af_reference_percent", "af_query_percent"]:
        frame[column] = pd.to_numeric(frame[column], errors="raise")
    if (frame["ani"] < 0).any() or (frame["ani"] > 100).any():
        raise ValueError("skani ANI values must lie in [0, 100]")
    for column in ["af_reference_percent", "af_query_percent"]:
        if (frame[column] < 0).any() or (frame[column] > 100).any():
            raise ValueError(f"skani {column} values must lie in [0, 100]")
    frame["af_reference"] = frame["af_reference_percent"] / 100.0
    frame["af_query"] = frame["af_query_percent"] / 100.0
    frame["af_max"] = frame[["af_reference", "af_query"]].max(axis=1)
    frame["reference_accession"] = frame["reference_fasta"].map(accession_from_path)
    if frame["reference_accession"].eq("").any():
        raise ValueError("Could not parse a GTDB reference accession from one or more Ref_file values")
    frame["query_path_norm"] = frame["query_fasta"].map(normalize_path)
    unknown_queries = sorted(set(frame["query_path_norm"]).difference(known_query_paths))
    if unknown_queries:
        raise ValueError(f"skani results contain a query absent from the query map: {unknown_queries[0]}")
    return frame


def load_metadata(path: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    metadata = read_table(path)
    required = {"accession", "gtdb_taxonomy", "ncbi_genome_category"}
    require_columns(metadata, required, "GTDB metadata")
    metadata = metadata.loc[:, list(required)].copy()
    metadata["accession_norm"] = metadata["accession"].map(normalize_accession)
    metadata["species"] = metadata["gtdb_taxonomy"].map(species_from_taxonomy)
    metadata["is_isolate"] = metadata["ncbi_genome_category"].str.strip().str.lower().eq("none")
    if metadata["accession_norm"].eq("").any() or metadata["species"].eq("").any():
        raise ValueError("GTDB metadata contains an empty normalized accession or species assignment")

    accession_species_counts = metadata.groupby("accession_norm")["species"].nunique()
    if (accession_species_counts > 1).any():
        raise ValueError("GTDB metadata maps an accession to multiple species")
    reference_metadata = metadata[["accession_norm", "species"]].drop_duplicates("accession_norm")
    species_summary = metadata.groupby("species", as_index=False).agg(
        species_genome_count=("accession_norm", "size"),
        species_isolate_count=("is_isolate", "sum"),
    )
    species_summary["species_has_isolate"] = species_summary["species_isolate_count"].gt(0)
    return reference_metadata, species_summary


def call_status(
    query_map: pd.DataFrame,
    skani: pd.DataFrame,
    reference_metadata: pd.DataFrame,
    species_summary: pd.DataFrame,
    ani_threshold: float,
    af_threshold: float,
) -> pd.DataFrame:
    annotated = skani.merge(
        reference_metadata,
        left_on="reference_accession",
        right_on="accession_norm",
        how="left",
        validate="many_to_one",
    )
    if annotated["species"].isna().any():
        accession = annotated.loc[annotated["species"].isna(), "reference_accession"].iloc[0]
        raise ValueError(f"GTDB metadata is missing reference accession: {accession}")
    annotated = annotated.merge(species_summary, on="species", how="left", validate="many_to_one")
    qualifying = annotated.loc[
        annotated["ani"].ge(ani_threshold) & annotated["af_max"].ge(af_threshold)
    ].copy()
    qualifying = qualifying.sort_values(
        ["query_path_norm", "species_has_isolate", "ani", "af_max", "reference_accession"],
        ascending=[True, False, False, False, True],
    )

    calls_by_query: dict[str, dict[str, object]] = {}
    for query_path in query_map["query_path_norm"].drop_duplicates():
        hits = qualifying.loc[qualifying["query_path_norm"].eq(query_path)]
        cultured_hits = hits.loc[hits["species_has_isolate"]]
        if not cultured_hits.empty:
            chosen = cultured_hits.iloc[0]
            status = "Cultured"
            reason = "qualifying species match contains at least one isolate genome"
        elif not hits.empty:
            chosen = hits.iloc[0]
            status = "Uncultured"
            reason = "qualifying species matches contain no isolate genome"
        else:
            chosen = None
            status = "Uncultured"
            reason = "no qualifying representative hit"

        calls_by_query[query_path] = {
            "status": status,
            "culture_call_reason": reason,
            "match_reference": str(chosen["reference_accession"]) if chosen is not None else "",
            "match_species": str(chosen["species"]) if chosen is not None else "",
            "match_ani": float(chosen["ani"]) if chosen is not None else "",
            "match_af_reference": float(chosen["af_reference"]) if chosen is not None else "",
            "match_af_query": float(chosen["af_query"]) if chosen is not None else "",
            "match_af_max": float(chosen["af_max"]) if chosen is not None else "",
            "species_isolate_count": int(chosen["species_isolate_count"]) if chosen is not None else 0,
            "qualifying_reference_hit_count": int(len(hits)),
            "qualifying_species_count": int(hits["species"].nunique()),
        }

    call_frame = pd.DataFrame.from_dict(calls_by_query, orient="index").rename_axis("query_path_norm").reset_index()
    output = query_map.merge(call_frame, on="query_path_norm", how="left", validate="many_to_one")
    output = output.drop(columns=["query_fasta", "query_path_norm"])
    if output["status"].isna().any() or not output["status"].isin(["Cultured", "Uncultured"]).all():
        raise ValueError("Culture-status calling produced a missing or invalid status")
    if not output.loc[output["status"].eq("Cultured"), "species_isolate_count"].gt(0).all():
        raise ValueError("A Cultured call lacks isolate evidence")
    return output


def write_itol(path: Path, calls: pd.DataFrame, cultured_color: str, uncultured_color: str) -> None:
    colors = {"Cultured": cultured_color, "Uncultured": uncultured_color}
    lines = [
        "DATASET_COLORSTRIP",
        "SEPARATOR\tTAB",
        "DATASET_LABEL\tCulture_status",
        "COLOR\t#000000",
        "COLOR_BRANCHES\t0",
        "LEGEND_TITLE\tStatus",
        "LEGEND_SHAPES\t1\t1",
        f"LEGEND_COLORS\t{cultured_color}\t{uncultured_color}",
        "LEGEND_LABELS\tCultured\tUncultured",
        "DATA",
    ]
    lines.extend(f"{row.tip_id}\t{colors[row.status]}" for row in calls.itertuples(index=False))
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--query-map", type=Path, required=True, help="TSV/CSV with tip_id and query_fasta")
    parser.add_argument("--skani-results", type=Path, required=True, help="Raw skani search table")
    parser.add_argument("--metadata", type=Path, required=True, help="GTDB metadata TSV or TSV.GZ")
    parser.add_argument("--status-output", type=Path, required=True, help="Output tip-level status/evidence TSV")
    parser.add_argument("--itol-output", type=Path, required=True, help="Output iTOL DATASET_COLORSTRIP file")
    parser.add_argument("--ani-threshold", type=float, default=95.0, help="Minimum ANI percentage")
    parser.add_argument("--af-threshold", type=float, default=0.30, help="Minimum max bidirectional AF fraction")
    parser.add_argument("--cultured-color", default=DEFAULT_CULTURED_COLOR)
    parser.add_argument("--uncultured-color", default=DEFAULT_UNCULTURED_COLOR)
    args = parser.parse_args()

    if not 0 <= args.ani_threshold <= 100:
        raise ValueError("--ani-threshold must lie in [0, 100]")
    if not 0 <= args.af_threshold <= 1:
        raise ValueError("--af-threshold must lie in [0, 1]")
    cultured_color = validate_color(args.cultured_color, "--cultured-color")
    uncultured_color = validate_color(args.uncultured_color, "--uncultured-color")

    query_map = load_query_map(args.query_map)
    skani = load_skani_results(args.skani_results, set(query_map["query_path_norm"]))
    reference_metadata, species_summary = load_metadata(args.metadata)
    calls = call_status(
        query_map, skani, reference_metadata, species_summary,
        args.ani_threshold, args.af_threshold,
    )
    args.status_output.parent.mkdir(parents=True, exist_ok=True)
    calls.to_csv(args.status_output, sep="\t", index=False)
    write_itol(args.itol_output, calls, cultured_color, uncultured_color)
    counts = calls["status"].value_counts().to_dict()
    print(f"Wrote {args.status_output} and {args.itol_output} (n={len(calls)}; {counts})")


if __name__ == "__main__":
    main()
