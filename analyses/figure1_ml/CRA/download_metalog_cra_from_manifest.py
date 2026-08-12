#!/usr/bin/env python3
"""Re-download available CRA Metalog files listed in the provenance manifest."""

from __future__ import annotations

import argparse
import csv
import gzip
import os
import subprocess
from datetime import datetime
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_MANIFEST = ROOT / "provenance" / "CRA_download_manifest.tsv"
OUTPUT_MANIFEST = ROOT / "CRA" / "download_run_manifest.tsv"
DOWNLOADABLE_STATUS = {"downloaded"}


def download(url: str, destination: Path, overwrite: bool) -> tuple[str, str]:
    destination.parent.mkdir(parents=True, exist_ok=True)
    if destination.exists() and not overwrite:
        return "already_exists", "not overwritten"
    temporary = destination.with_name(destination.name + ".part")
    if temporary.exists():
        temporary.unlink()
    command = [
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
        str(temporary),
    ]
    process = subprocess.run(command, capture_output=True, text=True)
    if process.returncode != 0:
        if temporary.exists():
            temporary.unlink()
        return "failed", process.stderr.strip() or f"curl exited {process.returncode}"
    if temporary.stat().st_size == 0:
        temporary.unlink()
        return "failed", "empty response"
    if destination.suffix == ".gz":
        try:
            with gzip.open(temporary, "rb") as handle:
                while handle.read(1024 * 1024):
                    pass
        except OSError as exc:
            temporary.unlink()
            return "failed", f"invalid gzip: {exc}"
    os.replace(temporary, destination)
    return "downloaded", "validated"


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--overwrite", action="store_true")
    args = parser.parse_args()

    rows: list[dict[str, object]] = []
    with args.manifest.open(newline="") as handle:
        for source in csv.DictReader(handle, delimiter="\t"):
            source_status = str(source.get("status", ""))
            destination = ROOT / str(source.get("local_file", ""))
            if source_status in DOWNLOADABLE_STATUS:
                status, message = download(str(source["url"]), destination, args.overwrite)
            else:
                status, message = "skipped", f"source manifest status: {source_status}"
            rows.append(
                {
                    **source,
                    "run_status": status,
                    "actual_bytes": destination.stat().st_size if destination.exists() else 0,
                    "message": message,
                    "checked_at": datetime.now().astimezone().isoformat(timespec="seconds"),
                }
            )

    OUTPUT_MANIFEST.parent.mkdir(parents=True, exist_ok=True)
    fields = list(rows[0]) if rows else ["study", "type", "url", "local_file", "status", "bytes"]
    with OUTPUT_MANIFEST.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote {OUTPUT_MANIFEST}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
