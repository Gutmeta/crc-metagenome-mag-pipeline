#!/usr/bin/env python3
"""Validate integrity and basic de-identification properties of the release."""

from __future__ import annotations

import csv
import gzip
import hashlib
import re
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "DATA_MANIFEST.tsv"
ALIAS_RE = re.compile(r"^SMP_[0-9a-f]{20}$")
METADATA_FIELDS = ["study_code", "sample_alias", "subject_disease_status", "age_years", "sex"]
PROFILE_FIELDS = {
    "metaphlan4_species": ["sample_alias", "clade_name", "rel_abund"],
}
FORBIDDEN_METADATA_FIELDS = {
    "subject_id",
    "longitude",
    "latitude",
    "collection_date",
    "sample_title",
    "sample_description",
    "medication",
}
FORBIDDEN_PREFIXES = ("/" + "mnt" + "/", "/" + "home" + "/")
TEXT_SUFFIXES = {".py", ".md", ".tsv", ".txt", ".sh"}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def validate_metadata(path: Path) -> tuple[int, int]:
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = reader.fieldnames or []
        if fields != METADATA_FIELDS:
            raise ValueError(f"{path.relative_to(ROOT)}: unexpected metadata fields {fields}")
        if FORBIDDEN_METADATA_FIELDS.intersection(fields):
            raise ValueError(f"{path.relative_to(ROOT)}: forbidden metadata field retained")
        aliases = []
        for row in reader:
            alias = row["sample_alias"]
            if not ALIAS_RE.fullmatch(alias):
                raise ValueError(f"{path.relative_to(ROOT)}: non-release alias {alias!r}")
            aliases.append(alias)
    if len(aliases) != len(set(aliases)):
        raise ValueError(f"{path.relative_to(ROOT)}: duplicate metadata aliases")
    return len(aliases), len(aliases)


def validate_profile(path: Path, file_type: str) -> tuple[int, int]:
    aliases: set[str] = set()
    rows = 0
    with gzip.open(path, "rt", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = reader.fieldnames or []
        if fields != PROFILE_FIELDS[file_type]:
            raise ValueError(f"{path.relative_to(ROOT)}: unexpected profile fields {fields}")
        for row in reader:
            alias = row["sample_alias"]
            if not ALIAS_RE.fullmatch(alias):
                raise ValueError(f"{path.relative_to(ROOT)}: non-release alias {alias!r}")
            float(row["rel_abund"])
            aliases.add(alias)
            rows += 1
    return len(aliases), rows


def scan_local_paths() -> list[str]:
    findings: list[str] = []
    for path in ROOT.rglob("*"):
        if not path.is_file() or path.suffix.lower() not in TEXT_SUFFIXES:
            continue
        try:
            text = path.read_text(errors="replace")
        except OSError:
            continue
        if any(prefix in text for prefix in FORBIDDEN_PREFIXES):
            findings.append(str(path.relative_to(ROOT)))
    return findings


def main() -> int:
    if not MANIFEST.exists():
        raise FileNotFoundError(MANIFEST)
    with MANIFEST.open(newline="") as handle:
        records = list(csv.DictReader(handle, delimiter="\t"))
    if not records:
        raise ValueError("DATA_MANIFEST.tsv contains no files")

    totals = {"files": 0, "samples": 0, "rows": 0, "bytes": 0}
    by_type: dict[str, int] = {}
    for record in records:
        relative = Path(record["relative_path"])
        if relative.is_absolute() or ".." in relative.parts:
            raise ValueError(f"Unsafe manifest path: {relative}")
        path = ROOT / relative
        if not path.is_file():
            raise FileNotFoundError(path)
        if path.stat().st_size != int(record["bytes"]):
            raise ValueError(f"{relative}: byte-count mismatch")
        if sha256(path) != record["sha256"]:
            raise ValueError(f"{relative}: SHA-256 mismatch")

        file_type = record["file_type"]
        if file_type == "metadata_minimized":
            sample_count, row_count = validate_metadata(path)
        elif file_type in PROFILE_FIELDS:
            sample_count, row_count = validate_profile(path, file_type)
        else:
            raise ValueError(f"Unknown manifest file type: {file_type}")
        if sample_count != int(record["samples"]) or row_count != int(record["rows"]):
            raise ValueError(f"{relative}: manifest sample/row-count mismatch")
        totals["files"] += 1
        totals["samples"] += sample_count
        totals["rows"] += row_count
        totals["bytes"] += path.stat().st_size
        by_type[file_type] = by_type.get(file_type, 0) + 1

    crosswalks = [path for path in ROOT.rglob("*") if path.is_file() and "crosswalk" in path.name.lower()]
    if crosswalks:
        raise ValueError("Identifier crosswalk-like files found: " + ", ".join(map(str, crosswalks)))
    path_findings = scan_local_paths()
    if path_findings:
        raise ValueError("Machine-specific absolute paths found in: " + ", ".join(path_findings))

    print("Public-release validation passed")
    print(f"Manifest files: {totals['files']} ({by_type})")
    print(f"Manifest rows: {totals['rows']}; bytes: {totals['bytes']}")
    print("Metadata schema, pseudonym format, file hashes, and path scan: PASS")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
