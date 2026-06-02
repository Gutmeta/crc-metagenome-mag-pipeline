#!/bin/bash

# 操作前改一下MIXED里面的前缀
REP_DIR=/path/to/crc-metagenome-mag-pipeline/C_genomes/drep_all_sa0.99_20260317_214714/final_reps
C_TCG_DIR=/path/to/crc-metagenome-mag-pipeline/C_genomes/C_TCG_genomes

mkdir -p "$C_TCG_DIR"

cp "$REP_DIR"/C1A/*.fa "$C_TCG_DIR"/
cp "$REP_DIR"/C1B/*.fa "$C_TCG_DIR"/

# 复制 MIXED，并把文件名前缀 C1A__ 或 C1B__ 改成 MIXED__
for f in "$REP_DIR"/MIXED/*.fa; do
  [ -e "$f" ] || continue   # 防止没有匹配文件时报错
  b=$(basename "$f")
  b=${b/#C1A__/MIXED__}
  b=${b/#C1B__/MIXED__}
  cp -- "$f" "$C_TCG_DIR/$b"
done