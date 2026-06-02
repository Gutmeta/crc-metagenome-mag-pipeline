#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import re
import glob
import string
import pandas as pd

# ============ 路径配置 ============
fa_dir = "/path/to/storage/data3/CRC/YachidaS_2019/Results/dRep_hq_bins_dir"
summary_path = "/path/to/crc-metagenome-mag-pipeline/classify_wf_out/YachidaS_2019/gtdbtk.summary.tsv"
group_path = "/path/to/crc-metagenome-mag-pipeline/Yachidas_2019_CRC_Group.txt"
output_path = "/path/to/storage/data3/CRC/YachidaS_2019/Results/gtdbtk_id_mapping_with_unclassified.tsv"

# 严格模式：True=发现规范化后重复则报错; False=自动为重复项添加 a/b/c… 后缀
STRICT_MODE = False

# ============ 工具函数 ============
def parse_sample_suffix_from_text(text: str):
    """
    从任意文本中解析 (Sample, genome_suffix)：
    匹配模式：([DES]RR数字+).(数字+)
    允许前缀(如 metabat2_)、后缀(_sub)和扩展名(.fa 等)。
    """
    m = re.search(r'([DES]RR\d+)\.(\d+)', str(text))
    if not m:
        return None, None
    return m.group(1), m.group(2)

def extract_best_taxonomy(class_str):
    """
    从 GTDB 分类串中，按 s__→g__→f__→o__→c__→p__→d__ 优先级提取第一个非空层级；
    若都缺失则返回 'Unclassified'
    """
    for rank in ['s__', 'g__', 'f__', 'o__', 'c__', 'p__', 'd__']:
        # 用 pandas 的提取以兼容 None/NaN
        match = pd.Series(class_str).str.extract(fr'{rank}([^;]*)').iloc[0, 0]
        if pd.notna(match) and match != "":
            return match
    return "Unclassified"

def disambiguate_suffixes(df):
    """
    自动消歧：对 (Sample, genome_suffix) 重复的行，按 fa_filename 排序后，
    将 genome_suffix 后追加 a/b/c…；返回消歧后的 DataFrame。
    """
    letters = list(string.ascii_lowercase)

    def _apply(group):
        if len(group) == 1:
            return group
        group = group.sort_values("fa_filename").copy()
        for i, idx in enumerate(group.index):
            tag = letters[i] if i < len(letters) else f"x{i}"
            group.at[idx, "genome_suffix"] = f"{group.at[idx, 'genome_suffix']}{tag}"
        return group

    return df.groupby(["Sample", "genome_suffix"], group_keys=False).apply(_apply)

# ============ 主流程 ============
def main():
    # 1) 读取 GTDB summary，并解析 Sample 与 genome_suffix（放宽匹配位置）
    gtdb_df = pd.read_csv(summary_path, sep="\t")
    gtdb_parse = gtdb_df["user_genome"].apply(lambda s: pd.Series(parse_sample_suffix_from_text(s)))
    gtdb_parse.columns = ["Sample", "genome_suffix"]
    gtdb_df = pd.concat([gtdb_df, gtdb_parse], axis=1)

    # 去除无法解析的行，避免后续合并产生歧义
    gtdb_df = gtdb_df.dropna(subset=["Sample", "genome_suffix"])
    gtdb_df["genome_suffix"] = gtdb_df["genome_suffix"].astype(str)

    # 提取 GTDB 最佳可用分类
    gtdb_df["GTDBTK"] = gtdb_df["classification"].apply(extract_best_taxonomy)

    # 2) 读取分组信息并准备 group_code 与 sample_suffix
    group_df = pd.read_csv(group_path, sep="\t", names=["Sample", "Group"], header=0)
    group_map = {"CRC": "CRC", "control": "CON"}
    group_df["group_code"] = group_df["Group"].map(group_map)
    group_df["sample_suffix"] = group_df["Sample"].astype(str).str[-3:]

    # 3) 扫描 fa 文件并解析出 Sample、genome_suffix
    fa_files = glob.glob(os.path.join(fa_dir, "*.fa"))
    records = []
    for path in fa_files:
        filename = os.path.basename(path)
        sample, suffix = parse_sample_suffix_from_text(filename)
        if not sample or not suffix:
            # 跳过无法解析的文件（也可改为记录日志）
            continue
        records.append({
            "Sample": sample,
            "genome_suffix": str(suffix),
            "fa_filename": filename,  # 保留原始文件名便于溯源
        })
    all_bins_df = pd.DataFrame(records)

    if all_bins_df.empty:
        raise RuntimeError("未在 fa_dir 中解析到任何合法的 .fa 文件，请检查命名或路径。")

    # 4) 规范化后检查 (Sample, genome_suffix) 是否重复
    dup_mask = all_bins_df.duplicated(subset=["Sample", "genome_suffix"], keep=False)
    dups = all_bins_df.loc[dup_mask].sort_values(["Sample", "genome_suffix", "fa_filename"])

    if not dups.empty:
        if STRICT_MODE:
            dup_report = os.path.join(os.path.dirname(output_path), "duplicate_bins_after_normalization.tsv")
            os.makedirs(os.path.dirname(dup_report), exist_ok=True)
            dups.to_csv(dup_report, sep="\t", index=False)
            raise ValueError(
                f"检测到规范化后的 (Sample, genome_suffix) 重复，共 {len(dups)} 条；"
                f"已写出清单：{dup_report}\n"
                f"样例：\n{dups.head(10)}"
            )
        else:
            all_bins_df = disambiguate_suffixes(all_bins_df)

    # 5) 合并分组与 GTDB 注释
    merged = all_bins_df.merge(
        group_df[["Sample", "group_code", "sample_suffix"]],
        on="Sample", how="left"
    )
    merged = merged.merge(
        gtdb_df[["Sample", "genome_suffix", "GTDBTK"]],
        on=["Sample", "genome_suffix"], how="left"
    )

    # 6) 填充未注释分类为 Unclassified，并构建统一 ID
    merged["GTDBTK"] = merged["GTDBTK"].fillna("Unclassified")
    merged["ID"] = merged["group_code"] + merged["sample_suffix"] + "." + merged["genome_suffix"]

    # 7) 输出
    final_df = merged[["ID", "GTDBTK", "fa_filename"]]
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    final_df.to_csv(output_path, sep="\t", index=False)

    # 终端提示
    print(f"[OK] 写出：{output_path}")
    if not dups.empty and not STRICT_MODE:
        print("[NOTE] 已对重复的 (Sample, genome_suffix) 自动消歧，请检查输出的一致性。")

if __name__ == "__main__":
    main()
