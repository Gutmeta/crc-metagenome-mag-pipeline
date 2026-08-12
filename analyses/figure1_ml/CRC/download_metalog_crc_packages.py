#!/usr/bin/env python3
"""Download and organize CRC Metalog packages by study code."""

from __future__ import annotations

import csv
import filecmp
import gzip
import json
import os
import shutil
import subprocess
import sys
from collections import Counter
from datetime import datetime
from pathlib import Path


BASE_URL = "https://metalog.embl.de"
PROFILE_TYPES = ("metaphlan4_species",)

ROOT = Path(__file__).resolve().parents[1]
CRC_DIR = ROOT / "CRC"
MANIFEST_PATH = CRC_DIR / "metalog_crc_download_manifest.tsv"
STUDIES_PATH = CRC_DIR / "metalog_crc_studies.tsv"

CRC_STUDIES = [
    {
        "study_code": "Feng_2015_CRC_Austria",
        "study_title": "Gut microbiome development along the colorectal adenoma-carcinoma sequence",
        "sample_count": 156,
        "legacy_dir": "Feng_2015",
    },
    {
        "study_code": "Gao_2021_CRC",
        "study_title": "Alterations, Interactions, and Diagnostic Potential of Gut Bacteria and Viruses in Colorectal Cancer",
        "sample_count": 126,
        "legacy_dir": "Gao_2021",
    },
    {
        "study_code": "Hannigan_2018_CRC",
        "study_title": "Diagnostic Potential and Interactive Dynamics of the Colorectal CRC Virome",
        "sample_count": 181,
        "legacy_dir": "Hannigan_2018",
    },
    {
        "study_code": "Liu_2022_CRC_China",
        "study_title": "Sequencing of genomic DNA extracted from stool samples of CRC and HS.",
        "sample_count": 165,
        "legacy_dir": "Liu_2022_CRC",
    },
    {
        "study_code": "Thomas_2019_CRC_Italy",
        "study_title": "Metagenomic analysis of colorectal cancer datasets identifies cross-cohort microbial diagnostic signatures and a link with choline degradation",
        "sample_count": 140,
        "legacy_dir": "Thomas_2019_CRC_Italy",
    },
    {
        "study_code": "Vogtmann_2016_CRC_USA",
        "study_title": "Reproducibility of associations between the human gut microbiome and colorectal cancer assessed in a patient population from Washington, DC, USA",
        "sample_count": 110,
        "legacy_dir": "Vogtmann_2016_CRC_USA",
    },
    {
        "study_code": "Wirbel_2019_CRC_Germany",
        "study_title": "Meta-analysis of fecal metagenomes reveals global microbial signatures that are specific for colorectal cancer",
        "sample_count": 120,
        "legacy_dir": "Wirbel_2019_CRC",
    },
    {
        "study_code": "Yachida_2019_CRC",
        "study_title": "Metagenomic and metabolomic analyses reveal distinct stage-specific phenotypes of the gut microbiota in colorectal cancer",
        "sample_count": 645,
        "legacy_dir": "Yachida_2019_CRC",
    },
    {
        "study_code": "Yu_2017_CRC_China",
        "study_title": "Metagenomic analysis of faecal microbiome as a tool towards targeted non-invasive biomarkers for colorectal cancer",
        "sample_count": 128,
        "legacy_dir": "Yu_2017",
    },
    {
        "study_code": "Zeller_2014_CRC_France",
        "study_title": "Potential of fecal microbiota for early-stage detection of colorectal cancer",
        "sample_count": 198,
        "legacy_dir": "Zeller_2014",
    },
]


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


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    with path.open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def metadata_counts(path: Path) -> tuple[str, bool]:
    counts: Counter[str] = Counter()
    with path.open("r", newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            status = (row.get("subject_disease_status") or "").strip()
            if status:
                counts[status] += 1
            else:
                counts["<blank>"] += 1
    contains_cra_group = any("adenoma" in label.lower() for label in counts)
    counts_text = "; ".join(f"{k}:{counts[k]}" for k in sorted(counts))
    return counts_text, contains_cra_group


def find_legacy_file(study_info: dict[str, object], package_type: str) -> Path | None:
    legacy_dir = CRC_DIR / str(study_info["legacy_dir"])
    if not legacy_dir.exists():
        return None

    if package_type == "metadata_all_wide":
        return None

    study = str(study_info["study_code"])
    patterns = [f"{package_type}_{study}_*.tsv.gz"]
    for pattern in patterns:
        matches = sorted(legacy_dir.glob(pattern))
        for candidate in matches:
            valid, _, _, _ = validate_file(candidate, package_type)
            if valid:
                return candidate
    return None


def files_identical(src: Path, dest: Path) -> bool:
    if not src.exists() or not dest.exists():
        return False
    try:
        return filecmp.cmp(src, dest, shallow=False)
    except OSError:
        return False


def copy_file(src: Path, dest: Path) -> str:
    dest.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(src, dest)
    return "existing_local_file"


def handle_package(study_info: dict[str, object], package_type: str) -> dict[str, object]:
    study = str(study_info["study_code"])
    study_dir = CRC_DIR / study
    dest = study_dir / package_filename(study, package_type)
    url = metadata_url(study) if package_type == "metadata_all_wide" else profile_url(study, package_type)
    legacy_src = find_legacy_file(study_info, package_type)

    rows = aliases = 0
    msg = "missing"

    if dest.exists():
        valid, msg, rows, aliases = validate_file(dest, package_type)
        if valid:
            if legacy_src and files_identical(legacy_src, dest):
                status = "existing_local_file"
            else:
                status = "already_exists"
        else:
            dest.unlink()
            status = "pending"
    else:
        status = "pending"

    if status == "pending":
        if legacy_src is not None:
            msg = copy_file(legacy_src, dest)
            status = "existing_local_file"
        else:
            ok, dl_msg = curl_file(url, dest)
            status = "downloaded" if ok else "unavailable"
            msg = dl_msg

    if status != "unavailable":
        valid, val_msg, rows, aliases = validate_file(dest, package_type)
        if not valid:
            status = "failed_validation"
            msg = val_msg

    return {
        "study_code": study,
        "study_title": study_info["study_title"],
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


def build_study_rows(manifest_rows: list[dict[str, object]]) -> list[dict[str, object]]:
    manifest_by_study: dict[str, list[dict[str, object]]] = {}
    for row in manifest_rows:
        manifest_by_study.setdefault(str(row["study_code"]), []).append(row)

    study_rows: list[dict[str, object]] = []
    for study_info in CRC_STUDIES:
        study = str(study_info["study_code"])
        study_dir = CRC_DIR / study
        metadata_path = study_dir / "metadata_all_wide.tsv"
        counts_text = ""
        contains_cra_group = False
        if metadata_path.exists():
            counts_text, contains_cra_group = metadata_counts(metadata_path)

        study_manifest = manifest_by_study.get(study, [])
        status_counter = Counter(str(row["status"]) for row in study_manifest)
        available_profiles = [
            str(row["package_type"])
            for row in study_manifest
            if row["package_type"] != "metadata_all_wide" and str(row["status"]) in {"already_exists", "existing_local_file", "downloaded"}
        ]
        unavailable_profiles = [
            str(row["package_type"])
            for row in study_manifest
            if row["package_type"] != "metadata_all_wide" and str(row["status"]) == "unavailable"
        ]
        study_rows.append(
            {
                "study_code": study,
                "study_title": study_info["study_title"],
                "sample_count": study_info["sample_count"],
                "legacy_dir": study_info["legacy_dir"],
                "standard_dir": str(study_dir),
                "contains_cra_group": "yes" if contains_cra_group else "no",
                "subject_disease_status_counts": counts_text,
                "available_profiles": ",".join(available_profiles),
                "unavailable_profiles": ",".join(unavailable_profiles),
                "status_counts": json.dumps(dict(sorted(status_counter.items())), ensure_ascii=True),
            }
        )
    return study_rows


def main() -> int:
    CRC_DIR.mkdir(parents=True, exist_ok=True)

    fields = [
        "study_code",
        "study_title",
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

    for study_info in CRC_STUDIES:
        study = str(study_info["study_code"])
        study_dir = CRC_DIR / study
        study_dir.mkdir(parents=True, exist_ok=True)
        source_rows: list[dict[str, object]] = []

        for package_type in ("metadata_all_wide", *PROFILE_TYPES):
            row = handle_package(study_info, package_type)
            source_rows.append(row)
            manifest_rows.append(row)

        write_tsv(study_dir / "source_info.tsv", source_rows, fields)

    write_tsv(MANIFEST_PATH, manifest_rows, fields)
    study_rows = build_study_rows(manifest_rows)
    write_tsv(
        STUDIES_PATH,
        study_rows,
        [
            "study_code",
            "study_title",
            "sample_count",
            "legacy_dir",
            "standard_dir",
            "contains_cra_group",
            "subject_disease_status_counts",
            "available_profiles",
            "unavailable_profiles",
            "status_counts",
        ],
    )

    print(f"CRC studies: {len(CRC_STUDIES)}")
    print(f"Manifest: {MANIFEST_PATH}")
    print("Status counts:", dict(Counter(str(row["status"]) for row in manifest_rows)))
    return 0


if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as exc:  # noqa: BLE001
        print(f"ERROR: {exc}", file=sys.stderr)
        raise
