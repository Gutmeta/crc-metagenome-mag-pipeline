#!/usr/bin/env python3
"""Download Metalog IBD-related study packages into dataset-named folders."""

from __future__ import annotations

import csv
import gzip
import json
import os
import re
import subprocess
import sys
from collections import Counter
from datetime import datetime
from pathlib import Path


BASE_URL = "https://metalog.embl.de"
STUDIES_URL = f"{BASE_URL}/dt_api/studies/all?draw=1&start=0&length=3000"
METADATA_TYPES = ("metadata_all_wide",)
PROFILE_TYPES = ("metaphlan4_species",)

ROOT = Path(__file__).resolve().parents[1]
IBD_DIR = ROOT / "IBD"
MANIFEST_PATH = IBD_DIR / "metalog_ibd_download_manifest.tsv"
STUDIES_PATH = IBD_DIR / "metalog_ibd_studies.tsv"


def now() -> str:
    return datetime.now().isoformat(timespec="seconds")


def curl_bytes(url: str) -> bytes:
    cmd = [
        "curl",
        "-4",
        "--retry",
        "5",
        "--retry-delay",
        "2",
        "--fail",
        "--location",
        "--show-error",
        "--silent",
        url,
    ]
    proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if proc.returncode != 0:
        raise RuntimeError(proc.stderr.decode("utf-8", "replace").strip() or f"curl exited {proc.returncode}")
    return proc.stdout


def curl_file(url: str, path: Path) -> tuple[bool, str]:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + ".part")
    if tmp.exists():
        tmp.unlink()
    cmd = [
        "curl",
        "-4",
        "--retry",
        "5",
        "--retry-delay",
        "2",
        "--fail",
        "--location",
        "--show-error",
        "--silent",
        url,
        "-o",
        str(tmp),
    ]
    proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if proc.returncode != 0:
        if tmp.exists():
            tmp.unlink()
        return False, proc.stderr.strip() or f"curl exited {proc.returncode}"
    os.replace(tmp, path)
    return True, "downloaded"


def is_ibd_study(row: dict[str, object]) -> bool:
    code = str(row.get("study_code") or "")
    title = str(row.get("study_title") or "")
    text = f"{code} {title}".lower()

    include = (
        "ibd" in text
        or "crohn" in text
        or "ulcerative colitis" in text
        or "inflammatory bowel" in text
        or re.search(r"(^|[_-])uc([_-]|$)", code.lower()) is not None
    )
    exclude = any(
        term in text
        for term in (
            "necrotizing enterocolitis",
            "premature infant",
            "premature_infant",
            "infant_nec",
            "_nec_",
        )
    )
    return include and not exclude


def fetch_ibd_studies() -> list[dict[str, object]]:
    payload = json.loads(curl_bytes(STUDIES_URL).decode("utf-8"))
    studies = [row for row in payload["data"] if is_ibd_study(row)]
    return sorted(studies, key=lambda r: str(r.get("study_code") or ""))


def gzip_valid(path: Path) -> tuple[bool, str]:
    try:
        with gzip.open(path, "rb") as fh:
            while fh.read(1024 * 1024):
                pass
        return True, "gzip_ok"
    except OSError as exc:
        return False, str(exc)


def validate_tsv(path: Path, gzipped: bool) -> tuple[bool, str, int, int]:
    """Return valid, message, data row count, unique sample_alias count."""
    opener = gzip.open if gzipped else open
    aliases: set[str] = set()
    rows = 0
    try:
        with opener(path, "rt", newline="") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            if not reader.fieldnames:
                return False, "missing_header", 0, 0
            if "sample_alias" not in reader.fieldnames:
                return False, f"missing_sample_alias_header:{','.join(reader.fieldnames)}", 0, 0
            for row in reader:
                rows += 1
                sample_alias = row.get("sample_alias")
                if sample_alias:
                    aliases.add(sample_alias)
        if rows == 0:
            return False, "empty_table", rows, len(aliases)
        return True, "validated", rows, len(aliases)
    except Exception as exc:  # noqa: BLE001
        return False, f"parse_failed:{exc}", 0, 0


def validate_file(path: Path, package_type: str) -> tuple[bool, str, int, int]:
    if not path.exists() or path.stat().st_size == 0:
        return False, "missing_or_empty", 0, 0
    gzipped = path.suffix == ".gz"
    if gzipped:
        ok, msg = gzip_valid(path)
        if not ok:
            return False, f"gzip_failed:{msg}", 0, 0
    valid, msg, rows, aliases = validate_tsv(path, gzipped=gzipped)
    if not valid:
        return valid, msg, rows, aliases
    if package_type == "metaphlan4_species":
        with gzip.open(path, "rt", newline="") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            checked = 0
            species_rows = 0
            for row in reader:
                checked += 1
                if str(row.get("clade_name", "")).startswith("s__"):
                    species_rows += 1
                if checked >= 1000 and species_rows:
                    break
            if species_rows == 0:
                return False, "no_s__species_rows_found", rows, aliases
    return True, msg, rows, aliases


def metadata_url(study: str, metadata_type: str) -> tuple[str, str]:
    _, scope, shape = metadata_type.split("_", 2)
    url = f"{BASE_URL}/api/study/{study}/metadata/human/{scope}.{shape.replace('_', '.')}.tsv"
    filename = f"{metadata_type}.tsv"
    return url, filename


def profile_url(study: str, profile_type: str) -> tuple[str, str]:
    filename = f"{profile_type}_{study}_latest.tsv.gz"
    url = f"{BASE_URL}/static/download/profiles/by_study/{filename}"
    return url, filename


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def handle_package(study: str, study_dir: Path, package_type: str, url: str, filename: str) -> dict[str, object]:
    path = study_dir / filename
    if path.exists():
        valid, msg, rows, aliases = validate_file(path, package_type)
        if valid:
            status = "already_exists"
        else:
            path.unlink()
            ok, msg = curl_file(url, path)
            status = "downloaded" if ok else "unavailable"
            rows = aliases = 0
    else:
        ok, msg = curl_file(url, path)
        status = "downloaded" if ok else "unavailable"
        rows = aliases = 0

    if status != "unavailable":
        valid, msg, rows, aliases = validate_file(path, package_type)
        if not valid:
            status = "failed_validation"

    return {
        "study_code": study,
        "package_type": package_type,
        "url": url,
        "local_file": str(path.relative_to(ROOT)),
        "status": status,
        "bytes": path.stat().st_size if path.exists() else 0,
        "row_count": rows,
        "sample_alias_count": aliases,
        "message": msg,
        "downloaded_at": now(),
    }


def main() -> int:
    studies = fetch_ibd_studies()
    study_fields = ("study_code", "study_title", "sample_count")
    write_tsv(
        STUDIES_PATH,
        [
            {
                "study_code": row.get("study_code") or "",
                "study_title": row.get("study_title") or "",
                "sample_count": row.get("sample_count") or "",
            }
            for row in studies
        ],
        list(study_fields),
    )

    manifest_rows: list[dict[str, object]] = []
    fields = [
        "study_code",
        "package_type",
        "url",
        "local_file",
        "status",
        "bytes",
        "row_count",
        "sample_alias_count",
        "message",
        "downloaded_at",
    ]

    for study_row in studies:
        study = str(study_row.get("study_code") or "")
        if not study:
            continue
        study_dir = IBD_DIR / study
        study_dir.mkdir(parents=True, exist_ok=True)
        source_rows: list[dict[str, object]] = []

        for package_type in METADATA_TYPES:
            url, filename = metadata_url(study, package_type)
            row = handle_package(study, study_dir, package_type, url, filename)
            source_rows.append(row)
            manifest_rows.append(row)

        for package_type in PROFILE_TYPES:
            url, filename = profile_url(study, package_type)
            row = handle_package(study, study_dir, package_type, url, filename)
            source_rows.append(row)
            manifest_rows.append(row)

        write_tsv(study_dir / "source_info.tsv", source_rows, fields)

    write_tsv(MANIFEST_PATH, manifest_rows, fields)
    status_counts = Counter(str(row["status"]) for row in manifest_rows)
    print(f"IBD studies: {len(studies)}")
    print(f"Manifest: {MANIFEST_PATH}")
    print("Status counts:", dict(status_counts))
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:  # noqa: BLE001
        print(f"ERROR: {exc}", file=sys.stderr)
        raise
