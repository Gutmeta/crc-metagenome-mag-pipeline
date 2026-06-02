#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse
import pandas as pd

# 你可以改成更接近论文的淡色
CULTURED_COLOR   = "#b3c7ff"   # 淡蓝（Cultured）
UNCULTURED_COLOR = "#c9f2c9"   # 淡绿（Uncultured）
MISSING_COLOR    = "#BDBDBD"   # 灰色（Missing/Unknown）

def smart_read_tsv(path):
    if path.endswith(".gz"):
        return pd.read_csv(path, sep="\t", dtype=str, compression="gzip", low_memory=False)
    return pd.read_csv(path, sep="\t", dtype=str, low_memory=False)

def pick_col(cols, candidates):
    for c in candidates:
        if c in cols:
            return c
    raise SystemExit(f"Cannot find columns among {candidates}. Available: {list(cols)[:50]} ...")

def norm_acc(x: str) -> str:
    if x is None:
        return ""
    x = str(x).strip()
    return x.replace("RS_", "").replace("GB_", "")

def status_from_ncbi_genome_category(v: str) -> str:
    # 更保守：None/缺失 => Missing（不要硬判 Cultured）
    if v is None:
        return "Cultured"
    s = str(v).strip().lower()
    if s in ("", "none", "nan"):
        return "Cultured"
    return "Uncultured"

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--summary",
                    default="/path/to/crc-metagenome-mag-pipeline/classify_wf_out/CC_TCG/gtdbtk.summary.tsv",
                    help="GTDB-Tk classify_wf summary.tsv")
    ap.add_argument("--metadata",
                    default="/mnt/data8/gtdbtk_data/release226/bac120_metadata_r226.tsv.gz",
                    help="GTDB bac120_metadata_r226.tsv.gz")
    ap.add_argument("--out",
                    default="/path/to/crc-metagenome-mag-pipeline/gtdbtk_denovo/itol_cultured_colorstrip.txt",
                    help="Output iTOL DATASET_COLORSTRIP file")
    args = ap.parse_args()

    # ---- 1) build tip -> status ----
    summ = smart_read_tsv(args.summary)
    user_col = pick_col(summ.columns, ["user_genome", "genome", "genome_id"])
    ref_col  = pick_col(summ.columns, ["closest_genome_reference", "closest_placement_reference"])
    summ["ref"] = summ[ref_col].map(norm_acc)

    meta = smart_read_tsv(args.metadata)
    acc_col = None
    for c in ["accession", "assembly_accession", "ncbi_assembly_accession"]:
        if c in meta.columns:
            acc_col = c
            break
    if acc_col is None:
        guess = [c for c in meta.columns if "accession" in c.lower()]
        if not guess:
            raise SystemExit("Cannot find an accession column in metadata; please check metadata header.")
        acc_col = guess[0]

    cat_col = pick_col(meta.columns, ["ncbi_genome_category"])
    meta["acc"] = meta[acc_col].map(norm_acc)
    meta = meta[["acc", cat_col]].drop_duplicates()
    acc2cat = dict(zip(meta["acc"], meta[cat_col]))

    tip2status = {}
    for _, r in summ.iterrows():
        tip = r[user_col]
        ref = r["ref"]
        if not tip:
            continue
        if ref:
            tip2status[tip] = status_from_ncbi_genome_category(acc2cat.get(ref, None))
        else:
            tip2status[tip] = "Missing"

    # ---- 2) write iTOL color strip with branch coloring ----
    out = []
    out += [
        "DATASET_COLORSTRIP",
        "SEPARATOR\tTAB",
        "DATASET_LABEL\tCultured_status",
        "COLOR\t#000000",
        "COLOR_BRANCHES\t0",
        "LEGEND_TITLE\tStatus",
        "LEGEND_SHAPES\t1\t1\t1",
        f"LEGEND_COLORS\t{CULTURED_COLOR}\t{UNCULTURED_COLOR}\t{MISSING_COLOR}",
        "LEGEND_LABELS\tCultured\tUncultured\tMissing",
        "DATA",
    ]

    for tip, st in tip2status.items():
        if st == "Cultured":
            c = CULTURED_COLOR
        elif st == "Uncultured":
            c = UNCULTURED_COLOR
        else:
            c = MISSING_COLOR
        out.append(f"{tip}\t{c}")

    with open(args.out, "w", encoding="utf-8") as w:
        w.write("\n".join(out) + "\n")

    print(f"Wrote {args.out} (n={len(tip2status)})")

if __name__ == "__main__":
    main()
