cd gtdbtk_denovo
grep -oE 'C1[AB]__[^:(),;]+' tree.C1AB.nwk | sort -u > tips.txt

awk -F'__' 'BEGIN{OFS="\t"}
{
  g=$1;
  c=(g=="C1A")?"#1f77b4":"#ff7f0e";
  print $0,c
}' tips.txt > id_color.tsv

python - <<'PY'
hdr = [
"DATASET_COLORSTRIP",
"SEPARATOR\tTAB",
"DATASET_LABEL\tTCG_guild",
"COLOR\t#000000",
"LEGEND_TITLE\tTCG",
"LEGEND_SHAPES\t1\t1",
"LEGEND_COLORS\t#1f77b4\t#ff7f0e",
"LEGEND_LABELS\tC1A\tC1B",
"DATA",
]
with open("itol_tcg_colorstrip.txt","w") as w:
    w.write("\n".join(hdr) + "\n")
    w.write(open("id_color.tsv").read())
print("Wrote itol_tcg_colorstrip.txt")
PY