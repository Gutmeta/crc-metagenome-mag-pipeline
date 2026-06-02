#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import re
from collections import Counter, defaultdict
from pathlib import Path

PALETTE = [
    "#4e79a7",
    "#f28e2b",
    "#59a14f",
    "#e15759",
    "#76b7b2",
    "#b07aa1",
    "#edc948",
    "#9c755f",
    "#bab0ac",
]

MISSING_COLOR = "#BDBDBD"
ACCESSION_RE = re.compile(r"\b(?:[SED]RR\d+|SAM[DN]\d+|ERS\d+|SRS\d+)\b")


def tip_to_accession(tip: str) -> str:
    name = tip.split("__", 1)[1] if "__" in tip else tip
    name = name.replace("_sub", "")
    return name.split(".", 1)[0]


def source_label(path: Path) -> str:
    label = path.stem.replace("_finalused", "")
    return re.sub(r"_[12]$", "", label)


def collect_accession_sources(metadata_dir: Path) -> dict[str, set[str]]:
    acc_to_sources = defaultdict(set)
    for path in sorted(metadata_dir.glob("*.txt")):
        source = source_label(path)
        text = path.read_text(encoding="utf-8", errors="ignore")
        for acc in ACCESSION_RE.findall(text):
            acc_to_sources[acc].add(source)
    return acc_to_sources


def write_colorstrip(path: Path, rows: list[tuple[str, str]]) -> None:
    sources = sorted({source for _, source in rows if source})
    color_map = {source: PALETTE[i % len(PALETTE)] for i, source in enumerate(sources)}
    color_map[""] = MISSING_COLOR

    out = [
        "DATASET_COLORSTRIP",
        "SEPARATOR\tTAB",
        "DATASET_LABEL\tDataset_origin",
        "COLOR\t#000000",
    ]
    if sources:
        out.extend([
            "LEGEND_TITLE\tDataset_origin",
            "LEGEND_SHAPES\t" + "\t".join(["1"] * len(sources)),
            "LEGEND_COLORS\t" + "\t".join(color_map[source] for source in sources),
            "LEGEND_LABELS\t" + "\t".join(sources),
        ])
    out.append("DATA")

    for tip, source in rows:
        out.append(f"{tip}\t{color_map.get(source, MISSING_COLOR)}")

    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(out) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Generate an iTOL DATASET_COLORSTRIP ring for source dataset/cohort."
    )
    parser.add_argument(
        "--tips",
        default="/path/to/crc-metagenome-mag-pipeline/gtdbtk_denovo/tips.txt",
        help="tips.txt, one tree tip per line.",
    )
    parser.add_argument(
        "--metadata_dir",
        default="/path/to/crc-metagenome-mag-pipeline/Metadata_Information",
        help="Directory containing cohort metadata text files.",
    )
    parser.add_argument(
        "--out",
        default="/path/to/crc-metagenome-mag-pipeline/gtdbtk_denovo/itol_dataset_origin_colorstrip.txt",
        help="Output iTOL DATASET_COLORSTRIP file.",
    )
    args = parser.parse_args()

    tips_path = Path(args.tips)
    metadata_dir = Path(args.metadata_dir)
    out_path = Path(args.out)

    tips = [line.strip() for line in tips_path.read_text(encoding="utf-8").splitlines() if line.strip()]
    acc_to_sources = collect_accession_sources(metadata_dir)

    rows = []
    missing = []
    ambiguous = []
    for tip in tips:
        acc = tip_to_accession(tip)
        sources = sorted(acc_to_sources.get(acc, []))
        if not sources:
            rows.append((tip, ""))
            missing.append((tip, acc))
        elif len(sources) > 1:
            rows.append((tip, sources[0]))
            ambiguous.append((tip, acc, sources))
        else:
            rows.append((tip, sources[0]))

    write_colorstrip(out_path, rows)

    counts = Counter(source for _, source in rows if source)
    print(f"Wrote {out_path} (n={len(rows)})")
    for source, count in sorted(counts.items()):
        print(f"  {source}: {count}")
    if missing:
        print(f"[WARN] Missing source for {len(missing)} tips: {missing[:10]}")
    if ambiguous:
        print(f"[WARN] Ambiguous source for {len(ambiguous)} tips; used first: {ambiguous[:10]}")


if __name__ == "__main__":
    main()
