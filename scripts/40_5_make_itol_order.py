#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Generate a single iTOL DATASET_COLORSTRIP ring at ORDER level (o__) from
GTDB-Tk classify_wf summary, with:
  - high-contrast palette
  - optional collapsing low-frequency orders into "Other"

Output:
  <out_prefix>_order.txt

Examples:
  # 低频(<5)归为 Other
  python make_itol_order_ring.py --min_count 5

  # 只保留Top10，其余归为 Other
  python make_itol_order_ring.py --top_n 10

  # 两者一起用
  python make_itol_order_ring.py --min_count 5 --top_n 10
"""

import argparse
import colorsys
from collections import Counter
import pandas as pd

PALETTES = {
    "set2": ["#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3",
             "#a6d854", "#ffd92f", "#e5c494", "#b3b3b3"],
    "dark2": ["#1b9e77", "#d95f02", "#7570b3", "#e7298a",
              "#66a61e", "#e6ab02", "#a6761d", "#666666"],
    "tab10": ["#4e79a7", "#f28e2b", "#e15759", "#76b7b2", "#59a14f",
              "#edc948", "#b07aa1", "#ff9da7", "#9c755f", "#bab0ac"],
    "tab20": ["#4e79a7", "#a0cbe8", "#f28e2b", "#ffbe7d", "#e15759",
              "#ff9d9a", "#76b7b2", "#59a14f", "#8cd17d", "#b6992d",
              "#f1ce63", "#499894", "#86bcb6", "#afafaf", "#d37295",
              "#fabfd2", "#b07aa1", "#d4a6c8", "#9d7660", "#d7b5a6"],
}

MISSING_COLOR = "#BDBDBD"   # taxonomy缺失/空
OTHER_COLOR   = "#808080"   # 低频合并后的Other（比缺失更深一点，便于区分）


def pick_col(cols, candidates):
    for c in candidates:
        if c in cols:
            return c
    raise SystemExit(f"Cannot find columns among {candidates}. Available: {list(cols)}")


def get_rank(tax, prefix):
    """Extract rank value from GTDB taxonomy string like d__;p__;c__;o__;f__;g__;s__"""
    if pd.isna(tax) or not tax:
        return ""
    for part in str(tax).split(";"):
        part = part.strip()
        if part.startswith(prefix):
            return part.replace(prefix, "")
    return ""


def golden_hsv_color(i, sat=0.70, val=0.90):
    """Extra distinct-ish colors using golden-ratio hue stepping."""
    h = (i * 0.61803398875) % 1.0
    r, g, b = colorsys.hsv_to_rgb(h, sat, val)
    return "#{:02x}{:02x}{:02x}".format(int(r * 255), int(g * 255), int(b * 255))


def make_color_map(categories, palette_name="tab20", missing=MISSING_COLOR):
    cats = sorted(set([c for c in categories if c and c != "Other"]))
    base = PALETTES.get(palette_name, PALETTES["tab20"])

    cmap = {"": missing, "Other": OTHER_COLOR}
    for idx, c in enumerate(cats):
        if idx < len(base):
            cmap[c] = base[idx]
        else:
            cmap[c] = golden_hsv_color(idx - len(base))
    return cmap


def write_colorstrip(path, label, rows, palette_name="tab20", legend_max=30):
    cats = sorted({c for _, c in rows if c})
    cmap = make_color_map(cats, palette_name=palette_name)

    out = []
    out.append("DATASET_COLORSTRIP")
    out.append("SEPARATOR\tTAB")
    out.append(f"DATASET_LABEL\t{label}")
    out.append("COLOR\t#000000")

    # Legend (avoid huge legends)
    cats_for_legend = [c for c in cats if c]  # exclude empty
    if cats_for_legend and len(cats_for_legend) <= legend_max:
        out.append(f"LEGEND_TITLE\t{label}")
        out.append("LEGEND_SHAPES\t" + "\t".join(["1"] * len(cats_for_legend)))
        out.append("LEGEND_COLORS\t" + "\t".join([cmap.get(c, MISSING_COLOR) for c in cats_for_legend]))
        out.append("LEGEND_LABELS\t" + "\t".join(cats_for_legend))

    out.append("DATA")
    for tip, cat in rows:
        out.append(f"{tip}\t{cmap.get(cat, MISSING_COLOR)}")

    with open(path, "w", encoding="utf-8") as w:
        w.write("\n".join(out) + "\n")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tips", default="/path/to/crc-metagenome-mag-pipeline/gtdbtk_denovo/tips.txt",
                    help="tips.txt (one tip per line, must match tree labels)")
    ap.add_argument("--summary", default="/path/to/crc-metagenome-mag-pipeline/classify_wf_out/CC_TCG/gtdbtk.summary.tsv",
                    help="GTDB-Tk classify_wf summary tsv")
    ap.add_argument("--out_prefix", default="/path/to/crc-metagenome-mag-pipeline/gtdbtk_denovo/itol",
                    help="output prefix (output will be <out_prefix>_order.txt)")

    ap.add_argument("--legend_max", type=int, default=30,
                    help="Max categories to show in legend")
    ap.add_argument("--order_palette", default="tab20", choices=list(PALETTES.keys()),
                    help="Palette for order ring (o__)")

    # collapsing options
    ap.add_argument("--min_count", type=int, default=5,
                    help="Collapse orders with count < min_count into 'Other' (0 = disable)")
    ap.add_argument("--top_n", type=int, default=10,
                    help="Keep only top N most frequent orders; rest -> 'Other' (0 = disable)")

    args = ap.parse_args()

    tips = [x.strip() for x in open(args.tips, encoding="utf-8") if x.strip()]

    df = pd.read_csv(args.summary, sep="\t", dtype=str)
    id_col = pick_col(df.columns, ["user_genome", "genome", "genome_id"])
    tax_col = pick_col(df.columns, ["classification", "gtdb_taxonomy", "taxonomy"])

    # tip -> taxonomy
    m = dict(zip(df[id_col], df[tax_col]))

    unmatched = [t for t in tips if t not in m]
    if unmatched:
        print(f"[WARN] {len(unmatched)} tips not found in summary (show 10): {unmatched[:10]}")

    # extract raw orders for all tips
    raw_orders = []
    for t in tips:
        tax = m.get(t, "")
        raw_orders.append(get_rank(tax, "o__"))

    # count non-empty
    cnt = Counter([o for o in raw_orders if o])
    kept = set(cnt.keys())

    # apply min_count
    if args.min_count and args.min_count > 0:
        kept = {o for o, c in cnt.items() if c >= args.min_count}

    # apply top_n among kept
    if args.top_n and args.top_n > 0:
        # sort by count desc
        ordered = sorted([(o, cnt[o]) for o in kept], key=lambda x: (-x[1], x[0]))
        kept = {o for o, _ in ordered[:args.top_n]}

    # build rows with collapsing
    rows = []
    other_n = 0
    for tip, o in zip(tips, raw_orders):
        if not o:
            rows.append((tip, ""))  # missing
        elif o in kept:
            rows.append((tip, o))
        else:
            rows.append((tip, "Other"))
            other_n += 1

    out_path = f"{args.out_prefix}_order.txt"
    write_colorstrip(
        path=out_path,
        label="order",
        rows=rows,
        palette_name=args.order_palette,
        legend_max=args.legend_max,
    )

    print("Wrote:")
    print(f"  - {out_path}")
    if args.min_count or args.top_n:
        print(f"[INFO] Collapsing enabled. Other tips: {other_n} / {len(tips)}")
        print(f"[INFO] Kept orders: {len(kept)} (missing not counted)")


if __name__ == "__main__":
    main()
