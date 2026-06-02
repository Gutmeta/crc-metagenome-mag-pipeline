#!/usr/bin/env bash
set -Eeuo pipefail

# FengQ_2015 的 204 本地结果保存在:
# /path/to/scratch/tmp_slurm/FengQ_2015/<sample>/<job_id>/output/dRep/dereplicated_genomes
# 同一样本可能有多次重跑，因此这里只处理“最新一次成功(job rc=0)”的结果目录。

base_dir="/path/to/scratch/tmp_slurm/FengQ_2015"
renamed=0
samples_checked=0

for sample_dir in "${base_dir}"/ERR*; do
    [ -d "$sample_dir" ] || continue
    sample=$(basename "$sample_dir")
    latest_job=""

    for meta in "$sample_dir"/*/meta.txt; do
        [ -f "$meta" ] || continue
        rc=$(awk -F= '/^rc=/{print $2}' "$meta")
        [ "$rc" = "0" ] || continue
        job=$(awk -F= '/^job_id=/{print $2}' "$meta")
        if [ -z "$latest_job" ] || [ "$job" -gt "$latest_job" ]; then
            latest_job="$job"
        fi
    done

    [ -n "$latest_job" ] || continue
    samples_checked=$((samples_checked + 1))

    genome_dir="$sample_dir/$latest_job/output/dRep/dereplicated_genomes"
    [ -d "$genome_dir" ] || continue

    for f in "$genome_dir"/*.fa; do
        [ -e "$f" ] || continue
        fname=$(basename "$f")
        if [[ "$fname" == bin*.fa ]]; then
            newname="${sample}${fname#bin}"
            mv "$f" "$genome_dir/$newname"
            echo "renamed: $f -> $genome_dir/$newname"
            renamed=$((renamed + 1))
        fi
    done
done

echo "checked_samples=${samples_checked}"
echo "renamed_files=${renamed}"
