# Claim Support Audit

| Claim | Decision | Evidence |
|---|---|---|
| CRA species-level classifier is lower than CRC and IBD | Supported | Exact permutation q values for CRA vs CRC/IBD are significant; see auc_disease_pairwise_tests.tsv. |
| CRC and IBD AUCs are high | Not supported | AUCs are moderate, with cohort heterogeneity. |
| Microbial AUC is similar to clinical indicators | Not assessed | No direct clinical-indicator model or biomarker AUC table was found in this analysis folder. |
| Four diseases share many associated species | Not supported | Top300 four-way shared_count is 0. |
| CRC and IBD associated species overlap | Supported | CRC|IBD Top300 q<0.05 and robust at larger thresholds. |
| IBS overlaps strongly with both CRC and IBD pairwise | Not supported as pairwise claim | CRC|IBS is not enriched; IBD|IBS is borderline in Top300. CRC|IBD|IBS three-way overlap is supported. |
