# Methods

## Study design and data acquisition

Human metadata and precomputed taxonomic profiles were obtained from the Metalog database (`https://metalog.embl.de/`). The principal retrieval and validation run was performed on 24 May 2026 (UTC+08:00); earlier locally cached profiles are identified in the accompanying provenance manifests. The primary analysis used Metalog MetaPhlAn 4 species-level relative-abundance tables (`metaphlan4_species_*_latest.tsv.gz`). Metadata and abundance records were linked using `sample_alias`. Downloaded files were checked for non-zero size, readable tabular or gzip structure, expected identifier columns, and the presence of species-level records.

## Disease definitions and cohort eligibility

Analyses were organized as disease-specific case-control arms. Controls were samples labelled `CTR` or `control patient`. CRA cases were samples whose disease-status label contained `adenoma` or `adenomatous`; CRC cases were samples labelled `colorectal cancer`, with adenoma samples excluded from the CRC comparison. IBD cases comprised Crohn's disease, ulcerative colitis, inflammatory bowel disease, or indeterminate colitis. IBS cases were labelled `irritable bowel syndrome`. Duplicate `sample_alias` records were resolved by retaining the first eligible record.

A disease arm was eligible when both metadata and a readable MetaPhlAn 4 species profile were available, sample identifiers overlapped, and at least 10 cases and 10 controls remained after matching. Duplicate CRA representations of the same study were resolved in favor of the standardized CRC-directory copy when available. The final arm-level dataset contained six CRA cohorts (615 cases and 1,102 controls), nine CRC cohorts (760 cases and 801 controls), nine IBD cohorts (2,271 cases and 791 controls), and one IBS cohort (323 cases and 151 controls). These are arm-level counts; several CRA and CRC arms from the same study reused the same control samples and were analyzed separately by disease.

## Species-abundance processing

Only MetaPhlAn 4 clades beginning with `s__` were retained. Long-format records containing `sample_alias`, `clade_name`, and `rel_abund` were pivoted to sample-by-species matrices; duplicate sample-species records were summed and absent observations were set to zero. Before association testing or model feature ranking, a species was retained if it was detected in at least 5% of samples or at least 10 samples, and its mean relative abundance was at least `1e-5`. Log-transformed analyses used `log10(relative abundance + 1e-6)`.

## Within-cohort association testing and machine-learning stability

Within each eligible disease arm, case and control abundances were compared for every retained species using a two-sided Mann-Whitney U test. Effect size was represented as `log2[(case mean + 1e-6)/(control mean + 1e-6)]`, and P values were adjusted within each arm by the Benjamini-Hochberg false-discovery-rate procedure. Differential support was defined as FDR q <= 0.10 and absolute log2 fold change >= 0.50.

Random-forest models were used to quantify complementary feature importance and stability. Stratified cross-validation used up to five folds, limited by the smaller class size. Each fold fitted 100 trees with `class_weight="balanced_subsample"`, square-root feature sampling, and a fixed random seed. Fold-level feature importance, rank, and frequency among the top 50 features were recorded. This random-forest analysis supported species ranking and was not used as the final cross-cohort classifier.

## Disease-level evidence synthesis and species ranking

Arm-level evidence was combined separately for each disease. For each species, two-sided P values were converted to signed Z scores using the sign of the arm-specific log2 fold change and combined by Stouffer's method with square-root sample-size weights. The disease-level log2 fold change was the correspondingly weighted mean. Meta-analysis P values were corrected within disease by Benjamini-Hochberg FDR. Directional consistency was the proportion of tested arms agreeing with the majority effect direction.

The ranking score was defined as `-log10(meta-q) * abs(meta-log2FC) * direction consistency * (1 + ML stability)`. Eligible ranked species required direction consistency >= 0.60 and evidence from at least two disease arms; one arm was sufficient for IBS, for which non-zero machine-learning stability was additionally required in the original Metalog analysis. The strict-species definition additionally required meta-q <= 0.10, absolute meta-log2FC >= 0.50, the required number of arm-level differential-support events, and high random-forest stability. Primary overlap analyses used broader ranked associated-species sets rather than the strict-species set.

Because only one primary IBS arm met the main inclusion criteria, the IBS ranking was augmented with standardized association results from the Goll_2020_FMT_IBS comparison and the published PRJEB34103 IBS-versus-healthy-control table. External source rows were eligible at source-level q <= 0.10 and absolute log2 fold change >= 0.50. Evidence across IBS sources was combined by the same signed-Z framework with square-root sample-size weights. Species names and effect directions were standardized so that positive effects consistently denoted enrichment in IBS.

## Species-panel classification and nested cohort validation

CRA, CRC, and IBD were modelled separately using class-weighted logistic regression (`class_weight="balanced"`, `C=1`, `liblinear`, maximum 1,000 iterations). Within each training set, species passing the abundance filter were ranked by `-log10(P) * abs(log2FC)`, where P was obtained from a two-sided Mann-Whitney U test. Ranked features were considered sequentially and a candidate was removed when its absolute Pearson correlation with any previously retained feature was >= 0.90. Candidate panel sizes ranged from 1 to 200.

Generalization performance was evaluated by nested leave-one-dataset-out (LODO) validation. In each outer fold, one complete cohort was withheld. Within the remaining cohorts, an inner LODO procedure selected the panel size with the highest mean validation ROC-AUC; ties were resolved in favor of the smaller panel. Feature filtering, ranking, correlation pruning, log transformation, and standardization were estimated from training data only. Standardization parameters were then applied unchanged to the held-out cohort. Performance was summarized by ROC-AUC. Cohort-specific 95% confidence intervals were obtained from 5,000 class-stratified bootstrap resamples, and one-sided label-permutation tests with 5,000 permutations assessed whether AUC exceeded 0.5.

The workflow also produced final fixed panels of 143 CRA, 18 CRC, and 185 IBD species. The final panel size was selected from the pooled inner-validation curves, and panel membership was derived by ranking the full disease-level dataset. Models were then refitted while leaving out each cohort. These fixed-panel AUCs are post-selection performance estimates because the feature identities were determined using all cohorts; strictly unbiased cross-cohort generalization claims therefore rely on the fully nested outer-LODO results.

## Comparison of predictive performance

Disease-level AUC distributions were summarized across held-out cohorts. CRA-versus-CRC and CRA-versus-IBD comparisons used one-sided exact permutation tests of the difference in mean cohort AUC, with the prespecified alternative that CRC or IBD performed better than CRA. CRC-versus-IBD used a two-sided exact permutation test. Mann-Whitney U tests were retained as supplementary comparisons. All disease-pair P values were adjusted by the Benjamini-Hochberg procedure. Shared-study paired sensitivity analyses were reported separately where applicable.

## Cross-disease overlap analysis

The principal association-overlap analysis selected the top 300 eligible species for CRA, CRC, and IBD and all 238 eligible standardized IBS species. For the final merged analysis, the CRA and CRC Top-300 sets were combined by union into a 533-species CRC/CRA set. Pairwise and multi-disease intersections were calculated among CRC/CRA, IBD, and IBS. These association sets were distinct from the smaller fixed prediction panels.

For each comparison, the null background was the union of species detectable in the participating diseases. Random sets matching the observed disease-specific set sizes were sampled without replacement from this background 10,000 times. The expected overlap was the mean null intersection size, and the empirical P value was `(1 + number of null overlaps >= observed overlap)/(10,000 + 1)`. Empirical P values were adjusted using Benjamini-Hochberg FDR. Enrichment was considered supported when the observed intersection exceeded the 95th percentile of the permutation null and permutation q < 0.05. Jaccard indices, overlap coefficients, observed-to-expected ratios, and pairwise hypergeometric tests were calculated as descriptive or supplementary statistics.

## Sensitivity and robustness analyses

Robustness was examined across Top-N thresholds of 50, 100, 150, 200, 300, and 500; after restricting to explicitly named species; after excluding individual studies shared between CRA and CRC; and after restricting CRC evidence to CRC-only studies. Where age and sex were sufficiently complete, species associations were additionally evaluated with binomial generalized linear models containing log-transformed abundance, age, and sex, followed by the same disease-level aggregation and permutation-overlap framework. Random seeds were fixed in all stochastic procedures.

## Software and reproducibility

Analyses were conducted in Python using NumPy, pandas, SciPy, scikit-learn, statsmodels, Matplotlib, and seaborn. Exact package versions are given in `requirements.txt` and `software_versions.tsv`. The `provenance/` directory records source URLs, original download or validation timestamps, availability status, file sizes, cohort inclusion decisions, and exclusions. A SHA-256 file manifest is supplied for integrity verification.
