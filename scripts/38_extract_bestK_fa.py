#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import shutil
import argparse
import glob

def read_filenames(txt_path: str):
    files = []
    with open(txt_path, "r", encoding="utf-8") as f:
        for line in f:
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            # 只取第一列（防止后面有分数/重要性等）
            fn = s.split()[0]
            files.append(fn)
    return files

def safe_basename(p: str) -> str:
    # 防止 txt 里带路径，统一只用 basename 做精确匹配
    return os.path.basename(p)

def main():
    ap = argparse.ArgumentParser(description="Strictly copy listed .fa files by exact filename.")
    ap.add_argument("--list",
                    default="/path/to/data2/CRC/CCDC2/ML_results/bestK_features_fullrange_xgb_gpu.txt",
                    help="Path to bestK_features_fullrange_xgb_gpu.txt")
    ap.add_argument("--src",
                    default="/path/to/crc-metagenome-mag-pipeline/C_genomes/C_TCG_genomes/",
                    help="Source dir containing fasta files")
    ap.add_argument("--dst",
                    default="/path/to/data2/CRC/CCDC2/CC_TCG_genomes",
                    help="Destination dir")
    ap.add_argument("--dry-run", action="store_true", help="Print actions only")
    args = ap.parse_args()

    if not os.path.isfile(args.list):
        print(f"[ERROR] list file not found: {args.list}", file=sys.stderr)
        sys.exit(1)
    if not os.path.isdir(args.src):
        print(f"[ERROR] source dir not found: {args.src}", file=sys.stderr)
        sys.exit(1)

    if args.dry_run:
        print(f"[DRY-RUN] would create dir: {args.dst}")
    else:
        os.makedirs(args.dst, exist_ok=True)

    filenames = read_filenames(args.list)

    copied = 0
    skipped_mixed = 0
    missing = []

    # 去重但保持顺序
    seen = set()
    uniq = []
    for fn in filenames:
        if fn not in seen:
            seen.add(fn)
            uniq.append(fn)

    for fn in uniq:
        base = safe_basename(fn)

        if base.upper().startswith("MIXED"):
            skipped_mixed += 1
            continue

        src_path = os.path.join(args.src, base)
        if not os.path.isfile(src_path):
            missing.append(base)
            continue

        dst_path = os.path.join(args.dst, base)
        if args.dry_run:
            print(f"[DRY-RUN] copy: {src_path} -> {dst_path}")
        else:
            shutil.copy2(src_path, dst_path)
        copied += 1

    # 删除目标目录里 MIXED*（不管后缀）
    mixed_files = glob.glob(os.path.join(args.dst, "MIXED*"))
    for p in mixed_files:
        if os.path.isfile(p):
            if args.dry_run:
                print(f"[DRY-RUN] delete: {p}")
            else:
                os.remove(p)

    miss_path = os.path.join(args.dst, "missing_files.txt")
    if args.dry_run:
        print(f"[DRY-RUN] would write missing list to: {miss_path} (n={len(missing)})")
    else:
        with open(miss_path, "w", encoding="utf-8") as f:
            for m in missing:
                f.write(m + "\n")

    print("====== DONE ======")
    print(f"List entries (unique): {len(uniq)}")
    print(f"Skipped MIXED*: {skipped_mixed}")
    print(f"Copied files: {copied}")
    print(f"Missing files: {len(missing)}")
    print(f"Missing list: {miss_path}")

if __name__ == "__main__":
    main()
