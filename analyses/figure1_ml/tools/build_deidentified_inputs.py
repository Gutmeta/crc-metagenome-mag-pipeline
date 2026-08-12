#!/usr/bin/env python3
"""Build minimized, pseudonymized analysis inputs from public Metalog exports.

No identifier crosswalk is written. Sample identifiers are replaced consistently
within a study by a SHA-256-derived release identifier. The output retains only
fields used by the released analysis: disease status, age, sex, and taxonomic
relative abundance.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import io
from pathlib import Path


DISEASE_ROOTS = ("CRA", "CRC", "IBD", "IBS")
METADATA_FIELDS = ("study_code", "sample_alias", "subject_disease_status", "age_years", "sex")
RELEASE_NAMESPACE = "Figure1_ML_public_release_v1"


def is_true(value: object) -> bool:
    return str(value).strip().lower() == "true"


def relative_input_path(value: str) -> Path:
    path = Path(value)
    if not path.is_absolute():
        if ".." in path.parts:
            raise ValueError(f"Unsafe relative input path {value!r}")
        return path
    for anchor in DISEASE_ROOTS:
        if anchor in path.parts:
            return Path(*path.parts[path.parts.index(anchor) :])
    raise ValueError(f"Cannot resolve public relative path from {value!r}")


def label_status(disease: str, status: str) -> str | None:
    normalized = str(status or "").strip().lower()
    if normalized in {"ctr", "control patient"}:
        return "control"
    if disease == "CRA":
        return "case" if ("adenoma" in normalized or "adenomatous" in normalized) else None
    if disease == "CRC":
        return "case" if normalized == "colorectal cancer" else None
    if disease == "IBD":
        return "case" if normalized in {
            "crohn's disease",
            "ulcerative colitis",
            "inflammatory bowel disease",
            "indeterminate colitis",
        } else None
    if disease == "IBS":
        return "case" if normalized == "irritable bowel syndrome" else None
    return None


def public_sample_id(study_key: str, source_identifier: str) -> str:
    payload = f"{RELEASE_NAMESPACE}|{study_key}|{source_identifier}".encode("utf-8")
    return "SMP_" + hashlib.sha256(payload).hexdigest()[:20]


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_included_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    return [row for row in rows if is_true(row.get("included_main", ""))]


def write_metadata(
    source_path: Path,
    destination: Path,
    diseases: set[str],
    study_key: str,
) -> tuple[dict[str, str], int]:
    destination.parent.mkdir(parents=True, exist_ok=True)
    identifier_map: dict[str, str] = {}
    seen: set[str] = set()
    written = 0
    with source_path.open(newline="") as source, destination.open("w", newline="") as target:
        reader = csv.DictReader(source, delimiter="\t")
        required = {"sample_alias", "subject_disease_status"}
        missing = required.difference(reader.fieldnames or [])
        if missing:
            raise ValueError(f"{source_path} missing columns {sorted(missing)}")
        writer = csv.DictWriter(target, fieldnames=METADATA_FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in reader:
            source_id = str(row.get("sample_alias", "")).strip()
            status = str(row.get("subject_disease_status", "")).strip()
            if not source_id or source_id in seen:
                continue
            if not any(label_status(disease, status) is not None for disease in diseases):
                continue
            seen.add(source_id)
            pseudonym = public_sample_id(study_key, source_id)
            identifier_map[source_id] = pseudonym
            writer.writerow(
                {
                    "study_code": row.get("study_code", destination.parent.name),
                    "sample_alias": pseudonym,
                    "subject_disease_status": status,
                    "age_years": row.get("age_years", ""),
                    "sex": row.get("sex", ""),
                }
            )
            written += 1
    if len(set(identifier_map.values())) != len(identifier_map):
        raise RuntimeError(f"Pseudonym collision in {study_key}")
    return identifier_map, written


def gzip_text_writer(path: Path):
    path.parent.mkdir(parents=True, exist_ok=True)
    raw = path.open("wb")
    compressed = gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0)
    text = io.TextIOWrapper(compressed, encoding="utf-8", newline="")
    return raw, compressed, text


def write_profile(
    source_path: Path,
    destination: Path,
    identifier_map: dict[str, str],
) -> tuple[int, int]:
    input_fields = ("sample_alias", "clade_name", "rel_abund")
    raw, compressed, target = gzip_text_writer(destination)
    written = 0
    samples: set[str] = set()
    try:
        with gzip.open(source_path, "rt", newline="") as source:
            reader = csv.DictReader(source, delimiter="\t")
            missing = set(input_fields).difference(reader.fieldnames or [])
            if missing:
                raise ValueError(f"{source_path} missing columns {sorted(missing)}")
            writer = csv.DictWriter(target, fieldnames=input_fields, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            for row in reader:
                source_id = str(row.get("sample_alias", ""))
                pseudonym = identifier_map.get(source_id)
                if pseudonym is None:
                    continue
                if not str(row.get("clade_name", "")).startswith("s__"):
                    continue
                output = {field: row.get(field, "") for field in input_fields}
                output["sample_alias"] = pseudonym
                writer.writerow(output)
                written += 1
                samples.add(pseudonym)
        target.flush()
    finally:
        target.close()
        if not compressed.closed:
            compressed.close()
        if not raw.closed:
            raw.close()
    return written, len(samples)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-root", type=Path, required=True)
    parser.add_argument("--output-root", type=Path, required=True)
    parser.add_argument("--main-inventory", type=Path, required=True)
    args = parser.parse_args()

    if args.source_root.resolve() == args.output_root.resolve():
        raise ValueError("Source and output roots must be different to protect source exports")

    inventories = [("metaphlan4_species", read_included_rows(args.main_inventory))]
    metadata_profiles: dict[Path, dict[str, object]] = {}
    for profile_type, rows in inventories:
        for row in rows:
            metadata_rel = relative_input_path(row["metadata_path"])
            profile_rel = relative_input_path(row["profile_path"])
            record = metadata_profiles.setdefault(metadata_rel, {"diseases": set(), "profiles": {}})
            record["diseases"].add(row["disease"])
            record["profiles"][profile_type] = profile_rel

    manifest_rows: list[dict[str, object]] = []
    for metadata_rel, record in sorted(metadata_profiles.items(), key=lambda item: str(item[0])):
        source_metadata = args.source_root / metadata_rel
        output_metadata = args.output_root / metadata_rel
        study_key = "/".join(metadata_rel.parts[:2])
        identifier_map, metadata_rows = write_metadata(
            source_metadata,
            output_metadata,
            set(record["diseases"]),
            study_key,
        )
        manifest_rows.append(
            {
                "disease_folder": metadata_rel.parts[0],
                "study_code": metadata_rel.parts[1],
                "file_type": "metadata_minimized",
                "relative_path": str(metadata_rel),
                "samples": metadata_rows,
                "rows": metadata_rows,
                "bytes": output_metadata.stat().st_size,
                "sha256": file_sha256(output_metadata),
                "retained_fields": ",".join(METADATA_FIELDS),
            }
        )
        for profile_type, profile_rel in sorted(record["profiles"].items()):
            source_profile = args.source_root / profile_rel
            output_profile = args.output_root / profile_rel
            rows_written, sample_count = write_profile(
                source_profile,
                output_profile,
                identifier_map,
            )
            manifest_rows.append(
                {
                    "disease_folder": profile_rel.parts[0],
                    "study_code": profile_rel.parts[1],
                    "file_type": profile_type,
                    "relative_path": str(profile_rel),
                    "samples": sample_count,
                    "rows": rows_written,
                    "bytes": output_profile.stat().st_size,
                    "sha256": file_sha256(output_profile),
                    "retained_fields": "sample_alias,clade_name,rel_abund",
                }
            )

    manifest_path = args.output_root / "DATA_MANIFEST.tsv"
    fields = [
        "disease_folder",
        "study_code",
        "file_type",
        "relative_path",
        "samples",
        "rows",
        "bytes",
        "sha256",
        "retained_fields",
    ]
    with manifest_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(manifest_rows)
    print(f"Wrote {len(manifest_rows)} files and {manifest_path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
