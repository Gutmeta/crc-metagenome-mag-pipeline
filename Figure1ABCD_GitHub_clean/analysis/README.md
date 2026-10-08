# Supporting calculations

`overlap/Overlap_results.tsv` reproduces Source Data Figure 1C. `overlap/Overlap_null_counts.npz` contains a 100000 × 4 integer matrix (`counts`), the corresponding comparison names (`comparisons`), `seed` and `iterations`. Every row represents one independent draw of all three sets from their respective backgrounds; columns contain CR–IBD, CR–IBS, IBD–IBS and the three-way intersection.

`auc_heterogeneity/AUC_descriptive_heterogeneity.tsv` reproduces the descriptive panel B source data. `AUC_exploratory_heterogeneity.tsv` contains logit-AUC and raw-AUC estimates of Q, I² and DerSimonian–Laird tau squared. `AUC_approximate_standard_errors.tsv` records the transformed values and CI-derived standard errors used in these calculations. The logit-AUC rows supply the I² labels in panel A.

See `docs/Figure1_methods.md` for formulas and interpretation limits, especially approximate variance recovery and dependence from overlapping LODO training sets.
