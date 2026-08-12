#!/usr/bin/env python3
"""Download and organize IBS Metalog packages by study code."""

from __future__ import annotations

import csv
import filecmp
import gzip
import os
import shutil
import subprocess
import sys
from collections import Counter
from datetime import datetime
from pathlib import Path


BASE_URL = "https://metalog.embl.de"
IBS_STUDIES = [
    {
        "study_code": "Goll_2020_FMT_IBS",
        "study_title": "Effects of fecal microbiota transplantation in subjects with irritable bowel syndrome are mirrored by changes in gut microbiome",
        "sample_count": 93,
    },
    {
        "study_code": "Mars_2020_IBS",
        "study_title": "Longitudinal Multi-omics Reveals Subset-Specific Mechanisms Underlying Irritable Bowel Syndrome",
        "sample_count": 481,
    },
    {
        "study_code": "PRJNA788943_IBS",
        "study_title": "IBS Metagenome",
        "sample_count": 56,
    },
]
PROFILE_TYPES = ("metaphlan4_species",)

ROOT = Path(__file__).resolve().parents[1]
IBS_DIR = ROOT / "IBS"
MANIFEST_PATH = IBS_DIR / "metalog_ibs_download_manifest.tsv"
STUDIES_PATH = IBS_DIR / "metalog_ibs_studies.tsv"

LOCAL_MARS_FILES = {
    "metadata_all_wide": IBS_DIR / "Mars_2020" / "metadata_Mars_2020_IBS_human_core_wide.tsv",
    "metaphlan4_species": IBS_DIR / "Mars_2020" / "metaphlan4_species_Mars_2020_IBS_2026-04-19.tsv.gz",
}


def now() -> str:
    return datetime.now().isoformat(timespec="seconds")


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


def gzip_valid(path: Path) -> tuple[bool, str]:
    try:
        with gzip.open(path, "rb") as fh:
            while fh.read(1024 * 1024):
                pass
        return True, "gzip_ok"
    except OSError as exc:
        return False, str(exc)


def validate_tsv(path: Path, gzipped: bool, package_type: str) -> tuple[bool, str, int, int]:
    opener = gzip.open if gzipped else open
    aliases: set[str] = set()
    rows = 0
    species_rows = 0
    try:
        with opener(path, "rt", newline="") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            if not reader.fieldnames:
                return False, "missing_header", 0, 0
            fields = set(reader.fieldnames)
            if "sample_alias" not in fields:
                return False, f"missing_sample_alias_header:{','.join(reader.fieldnames)}", 0, 0
            if package_type == "metadata_all_wide" and "subject_disease_status" not in fields:
                return False, "metadata_missing_subject_disease_status", 0, 0
            if package_type in {"metaphlan4", "metaphlan4_species"}:
                required = {"sample_alias", "clade_name", "rel_abund"}
                if not required.issubset(fields):
                    return False, f"bad_metaphlan_header:{','.join(reader.fieldnames)}", 0, 0
            for row in reader:
                rows += 1
                sample_alias = row.get("sample_alias")
                if sample_alias:
                    aliases.add(sample_alias)
                if package_type == "metaphlan4_species" and str(row.get("clade_name", "")).startswith("s__"):
                    species_rows += 1
        if rows == 0:
            return False, "empty_table", rows, len(aliases)
        if package_type == "metaphlan4_species" and species_rows == 0:
            return False, "no_s__species_rows_found", rows, len(aliases)
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
    return validate_tsv(path, gzipped, package_type)


def metadata_url(study: str) -> str:
    return f"{BASE_URL}/api/study/{study}/metadata/human/all.wide.tsv"


def profile_url(study: str, profile_type: str) -> str:
    return f"{BASE_URL}/static/download/profiles/by_study/{profile_type}_{study}_latest.tsv.gz"


def package_filename(study: str, package_type: str) -> str:
    if package_type == "metadata_all_wide":
        return "metadata_all_wide.tsv"
    return f"{package_type}_{study}_latest.tsv.gz"


def copy_local_if_valid(study: str, package_type: str, dest: Path) -> tuple[bool, str]:
    if study != "Mars_2020_IBS" or package_type not in LOCAL_MARS_FILES:
        return False, "no_local_cache"
    src = LOCAL_MARS_FILES[package_type]
    valid, msg, _, _ = validate_file(src, package_type)
    if not valid:
        return False, f"local_cache_invalid:{msg}"
    dest.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src, dest)
    return True, "existing_local_file"


def is_existing_local_file(study: str, package_type: str, dest: Path) -> bool:
    if study != "Mars_2020_IBS" or package_type not in LOCAL_MARS_FILES:
        return False
    src = LOCAL_MARS_FILES[package_type]
    if not src.exists() or not dest.exists():
        return False
    try:
        return filecmp.cmp(src, dest, shallow=False)
    except OSError:
        return False


def handle_package(study: str, package_type: str) -> dict[str, object]:
    study_dir = IBS_DIR / study
    dest = study_dir / package_filename(study, package_type)
    url = metadata_url(study) if package_type == "metadata_all_wide" else profile_url(study, package_type)

    if dest.exists():
        valid, msg, rows, aliases = validate_file(dest, package_type)
        if valid:
            status = "existing_local_file" if is_existing_local_file(study, package_type, dest) else "already_exists"
        else:
            dest.unlink()
            status = "pending"
    else:
        status = "pending"
        rows = aliases = 0
        msg = "missing"

    if status == "pending":
        copied, copy_msg = copy_local_if_valid(study, package_type, dest)
        if copied:
            status = "existing_local_file"
            msg = copy_msg
        else:
            ok, dl_msg = curl_file(url, dest)
            status = "downloaded" if ok else "unavailable"
            msg = dl_msg

    if status != "unavailable":
        valid, val_msg, rows, aliases = validate_file(dest, package_type)
        if not valid:
            status = "failed_validation"
            msg = val_msg
        elif msg in {"downloaded", "existing_local_file", "validated"}:
            msg = val_msg if status == "already_exists" else msg

    return {
        "study_code": study,
        "package_type": package_type,
        "url": url,
        "local_file": str(dest.relative_to(ROOT)),
        "status": status,
        "bytes": dest.stat().st_size if dest.exists() else 0,
        "row_count": rows,
        "sample_alias_count": aliases,
        "message": msg,
        "downloaded_at": now(),
    }


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main() -> int:
    IBS_DIR.mkdir(parents=True, exist_ok=True)
    write_tsv(STUDIES_PATH, IBS_STUDIES, ["study_code", "study_title", "sample_count"])

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
    manifest_rows: list[dict[str, object]] = []

    for study_info in IBS_STUDIES:
        study = study_info["study_code"]
        study_dir = IBS_DIR / study
        study_dir.mkdir(parents=True, exist_ok=True)
        source_rows: list[dict[str, object]] = []
        for package_type in ("metadata_all_wide", *PROFILE_TYPES):
            row = handle_package(study, package_type)
            source_rows.append(row)
            manifest_rows.append(row)
        write_tsv(study_dir / "source_info.tsv", source_rows, fields)

    write_tsv(MANIFEST_PATH, manifest_rows, fields)
    print(f"IBS studies: {len(IBS_STUDIES)}")
    print(f"Manifest: {MANIFEST_PATH}")
    print("Status counts:", dict(Counter(str(row["status"]) for row in manifest_rows)))
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:  # noqa: BLE001
        print(f"ERROR: {exc}", file=sys.stderr)
        raise
