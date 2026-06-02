#!/usr/bin/env python3

import argparse
import csv
from pathlib import Path

DATASETS = [
    "FengQ_2015",
    "ONCOBIOME_2025",
    "YachidaS_2019",
    "YangJ_2020",
    "YuJ_2015",
]


def taxa_from_ref_paths(ref_paths: Path) -> list[str]:
    taxa = []
    with ref_paths.open(encoding="utf-8") as handle:
        for line in handle:
            path = line.strip()
            if path:
                taxa.append(Path(path).name)
    if not taxa:
        raise SystemExit(f"No taxa found in ref paths file: {ref_paths}")
    return taxa


def selected_samples(group_file: Path, groups: set[str]) -> list[str]:
    samples = []
    seen = set()
    with group_file.open(encoding="utf-8", errors="ignore", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or "Sample" not in reader.fieldnames or "Group" not in reader.fieldnames:
            raise SystemExit(f"{group_file} must contain Sample and Group columns.")
        for row in reader:
            sample = (row.get("Sample") or "").strip()
            group = (row.get("Group") or "").strip().lower()
            if sample and group in groups and sample not in seen:
                samples.append(sample)
                seen.add(sample)
    return sorted(samples)


def read_abundance(path: Path) -> dict[str, str]:
    values = {}
    with path.open(encoding="utf-8", errors="ignore", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or "taxa.name" not in reader.fieldnames or "count.estimate" not in reader.fieldnames:
            raise SystemExit(f"{path} must contain taxa.name and count.estimate columns.")
        for row in reader:
            taxon = (row.get("taxa.name") or "").strip()
            value = (row.get("count.estimate") or "0").strip()
            if taxon:
                values[taxon] = value if value else "0"
    return values


def write_matrix(out_path: Path, taxa: list[str], samples: list[str], abundance_dir: Path) -> dict[str, int]:
    sample_values = {}
    missing_files = []
    for sample in samples:
        path = abundance_dir / f"abundance_{sample}.txt"
        if path.is_file():
            sample_values[sample] = read_abundance(path)
        else:
            missing_files.append(sample)

    out_samples = sorted(sample_values)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["taxa.name", *out_samples])
        for taxon in taxa:
            writer.writerow([taxon, *[sample_values[sample].get(taxon, "0") for sample in out_samples]])

    return {
        "selected_samples": len(samples),
        "written_samples": len(out_samples),
        "missing_files": len(missing_files),
        "taxa": len(taxa),
    }


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Build CRC/control sample C_TCG DiTASiC abundance matrices in wide TSV format."
    )
    parser.add_argument(
        "--ref_paths",
        default="/path/to/data2/CRC/CCDC2/DiTASiC/ref_paths.txt",
        help="C_TCG DiTASiC ref_paths.txt used for row order.",
    )
    parser.add_argument(
        "--data_root",
        default="/path/to/data2/CRC/CCDC2",
        help="Root containing <dataset>/DiTASiC/abundance.",
    )
    parser.add_argument(
        "--group_dir",
        default="/path/to/crc-metagenome-mag-pipeline",
        help="Directory containing <dataset>_CRC_Group.txt files.",
    )
    parser.add_argument(
        "--out_dir",
        default="/path/to/crc-metagenome-mag-pipeline/CRC_C_TCG_abundance_matrices",
        help="Output directory for *_CRC_C_TCG_abundance_matrix.tsv files.",
    )
    parser.add_argument(
        "--datasets",
        nargs="+",
        default=DATASETS,
        help="Dataset names to process.",
    )
    parser.add_argument(
        "--groups",
        nargs="+",
        default=["crc", "control"],
        help="Group labels to include, matched case-insensitively.",
    )
    args = parser.parse_args()

    taxa = taxa_from_ref_paths(Path(args.ref_paths))
    data_root = Path(args.data_root)
    group_dir = Path(args.group_dir)
    out_dir = Path(args.out_dir)
    groups = {g.lower() for g in args.groups}

    summary_rows = []
    for dataset in args.datasets:
        group_file = group_dir / f"{dataset}_CRC_Group.txt"
        abundance_dir = data_root / dataset / "DiTASiC" / "abundance"
        if not group_file.is_file():
            raise SystemExit(f"Missing group file: {group_file}")
        if not abundance_dir.is_dir():
            raise SystemExit(f"Missing abundance directory: {abundance_dir}")

        samples = selected_samples(group_file, groups)
        out_path = out_dir / f"{dataset}_CRC_C_TCG_abundance_matrix.tsv"
        stats = write_matrix(out_path, taxa, samples, abundance_dir)
        summary_rows.append((dataset, out_path, abundance_dir, stats))

    manifest_path = out_dir / "CRC_C_TCG_abundance_matrix_manifest.tsv"
    with manifest_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(
            [
                "dataset",
                "matrix_path",
                "source_abundance_dir",
                "taxa",
                "selected_samples",
                "written_samples",
                "missing_selected_abundance_files",
            ]
        )
        for dataset, out_path, abundance_dir, stats in summary_rows:
            writer.writerow(
                [
                    dataset,
                    out_path,
                    abundance_dir,
                    stats["taxa"],
                    stats["selected_samples"],
                    stats["written_samples"],
                    stats["missing_files"],
                ]
            )

    for dataset, out_path, _abundance_dir, stats in summary_rows:
        print(
            f"{dataset}: wrote {out_path} "
            f"({stats['taxa']} taxa x {stats['written_samples']} selected samples; "
            f"{stats['missing_files']} selected samples without abundance file)"
        )
    print(f"manifest: wrote {manifest_path}")


if __name__ == "__main__":
    main()
