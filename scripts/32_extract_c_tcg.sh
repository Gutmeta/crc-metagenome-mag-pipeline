#!/bin/bash

# Source and destination directories for dereplicated genomes.
REP_DIR=/path/to/drep_output/final_reps
C_TCG_DIR=/path/to/crc-metagenome-mag-pipeline/C_genomes/C_TCG_genomes

mkdir -p "$C_TCG_DIR"

cp "$REP_DIR"/C1A/*.fa "$C_TCG_DIR"/
cp "$REP_DIR"/C1B/*.fa "$C_TCG_DIR"/

# Copy MIXED genomes and normalize C1A__/C1B__ prefixes to MIXED__.
for f in "$REP_DIR"/MIXED/*.fa; do
  [ -e "$f" ] || continue   # Skip when the glob has no matches.
  b=$(basename "$f")
  b=${b/#C1A__/MIXED__}
  b=${b/#C1B__/MIXED__}
  cp -- "$f" "$C_TCG_DIR/$b"
done
