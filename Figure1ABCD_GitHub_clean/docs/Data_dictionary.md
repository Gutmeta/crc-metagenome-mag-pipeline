# Data dictionary

## `input/panel_a_lodo_auc.tsv`

- `disease`: analysis group (CRA, CRC or IBD).
- `project_accession`: held-out INSDC BioProject accession displayed in Figure 1A.
- `auc`: held-out-project ROC-AUC.
- `ci_lower`, `ci_upper`: lower and upper limits of the stratified bootstrap 95% confidence interval.
- `heldout_cases`, `heldout_controls`: numbers of cases and controls in the outer test project.
- `selected_k`: feature number selected by inner LODO validation within the outer training set.
- `bootstrap_iterations`: number of stratified bootstrap replicates.
- `permutation_p_one_sided`: empirical one-sided label-permutation *p* value for AUC at least as large as observed.
- `permutation_iterations`: number of label permutations.

## `input/species_sets.tsv`

- `analysis_set`: `CRC_CRA`, `IBD` or `IBS_Mars_2020`.
- `species`: species-level taxon identifier included in a selected species set. Membership is not an individual-species significance flag. The IBS set contains all 81 eligible candidates selected in the random-forest Top 50 in at least one of five folds.

## `input/species_detection_backgrounds.tsv`

- `analysis_set`: disease-analysis set associated with the detection background.
- `species`: species eligible to be sampled in fixed-size overlap permutations. The 465-species IBS background is the result of the abundance and detection filter described in Methods, before random-forest candidate selection.

## `input/ibs_species_selection.tsv`

One row for each of the 465 filtered Mars_2020 IBS background species. This table contains species-level selection evidence and aggregate statistics, not individual sample records. Abundances retain the source-export relative-abundance scale.

- `study_code`: source cohort code, `Mars_2020_IBS`.
- `species`: species-level taxon identifier, matching the IBS detection background.
- `case_samples`, `control_samples`: abundance-profile denominators, 323 and 151. These are sample counts, not confirmed counts of independent participants, and are repeated in every row.
- `detected_samples`: number of the 474 case and control samples with abundance greater than zero for the species.
- `detection_fraction`: `detected_samples / 474`.
- `mean_relative_abundance`: mean abundance across all 474 samples, including zero values for missing sample-species records.
- `case_mean_relative_abundance`, `control_mean_relative_abundance`: mean abundance within the respective sample group, including zeros.
- `log2_fold_change`: `log2((case_mean_relative_abundance + 10^-6) / (control_mean_relative_abundance + 10^-6))`; positive values indicate greater case mean abundance.
- `p_value`: two-sided Mann–Whitney p value comparing case and control abundances.
- `bh_q_value`: Benjamini–Hochberg-adjusted differential-test p value across the 465 IBS background species. It is distinct from the set-overlap `permutation_bh_q` in Source Data Figure 1C.
- `rf_top50_selection_frequency`: fraction of five stratified random-forest folds in which the species had importance rank at most 50 in the training model. Values range from 0 to 1 in increments of 0.2; this is a selection frequency, not a significance probability.
- `direction_consistency`: fraction of nonzero comparison-level fold-change directions matching the majority direction. For this single-cohort comparison, a value of 1 does not indicate external replication.
- `tested_comparisons`: number of eligible comparisons testing the species; 1 for this single-cohort IBS analysis.
- `eligible_for_selection`: whether direction consistency is at least 0.60, `tested_comparisons` is at least 1, and `rf_top50_selection_frequency` is greater than zero. No hard q-value or absolute fold-change threshold is included.
- `in_figure1_set`: membership in the 81-species IBS set used in Figure 1C–D. All 81 eligible candidates were retained because fewer than 300 were eligible.

## `Figure1_cohort_species_public.zip`

The archive contains sample-level profiles and minimal disease-comparison grouping files for 20 cohorts. Each `profiles/*.species.tsv.gz` contains `sample_alias` (pseudonymized release sample key), `clade_name` (species identifier beginning with `s__`) and `rel_abund` (source-export relative abundance). The key matches `sample_id` in the corresponding `sample_groups/*.groups.tsv`. Applicable `CRA`, `CRC`, `CRC_CRA`, `IBD` and `IBS` columns contain `case`, `control` or `not_in_comparison`.

Samples can include repeated visits or longitudinal specimens, and some controls are reused across comparisons. The archive README and cohort tables describe the source URLs, sample counts and file hashes. These processed species profiles are not sequencing reads or absolute organism counts.

## `source_data/Source_Data_Figure1A.tsv`

Contains the numeric values underlying Figure 1A, including exact AUCs, confidence intervals, sample counts, selected *K* and permutation settings.

## `source_data/Source_Data_Figure1B.tsv`

Contains the number of projects, minimum, first quartile, median, third quartile, maximum, interquartile range and total range for each disease. All cohorts have equal weight.

- `mean_auc`: arithmetic mean of cohort AUCs.
- `standard_deviation_auc`: sample standard deviation of cohort AUCs (n−1 denominator).
- `interquartile_range`: Q3−Q1, represented by the thick segments in Figure 1B.

Exploratory CI-derived I², Q and τ², on logit-AUC and raw-AUC scales, are kept in `analysis/auc_heterogeneity/`; their formulas and assumptions are documented there.

## `source_data/Source_Data_Figure1C.tsv`

- `comparison_key`: participating analysis sets joined by `|`, in the same order as size fields.
- `comparison`: displayed set comparison.
- `number_of_sets`, `set_sizes`: number and sizes of participating selected sets.
- `own_background_sizes`: corresponding group-specific detectable-background sizes, joined by `|`; there is no shared sampling pool.
- `common_background_species`: size of the intersection of participating detection backgrounds, used to calculate the analytical expectation (not the sampling pool).
- `observed_shared_species`: observed selected-set intersection size.
- `theoretical_expected`: exact expectation under independent uniform sampling from each group's own background.
- `permutation_expected`: mean null overlap from 100,000 own-background random draws.
- `expected_monte_carlo_se`: sample standard deviation of null overlaps divided by the square root of the number of draws.
- `permutation_95th_percentile`: 95th percentile of the simulated null.
- `observed_expected_ratio`: observed overlap divided by simulated expectation; displayed above Figure 1C bars as fold enrichment (×).
- `observed_theoretical_expected_ratio`: observed overlap divided by analytical expectation.
- `exceedance_count`: number of simulated overlaps greater than or equal to the observed overlap.
- `permutation_iterations`, `random_seed`: simulation count (100000) and seed (20260903).
- `permutation_p`: empirical upper-tail p value, `(exceedance_count + 1) / (permutation_iterations + 1)`.
- `tail_probability_mc_95ci_lower`, `tail_probability_mc_95ci_upper`: exact 95% binomial interval for the underlying null tail probability; Monte Carlo uncertainty only.
- `permutation_bh_q`: Benjamini–Hochberg-adjusted overlap p value across the four set comparisons; not an individual-species differential-test q value.
- `significant_enrichment`: whether q < 0.05 and observed overlap exceeds the null 95th percentile.

## `source_data/Source_Data_Figure1D.tsv`

- `row_set`, `column_set`: row and column identities in the symmetric matrix.
- `shared_species`: selected-set intersection size; diagonal values are set sizes.

## `source_data/Source_Data_Figure1A_heterogeneity.tsv`

Contains the exploratory logit-AUC heterogeneity estimates annotated in Figure 1A, also in the workbook sheet `Figure 1A heterogeneity`.

- `exploratory_i_squared_percent`: CI-derived I² percentage; rounded to one decimal place in the figure.
- `cochran_q`, `degrees_of_freedom`, `exploratory_q_p`: approximate Cochran Q and its degrees of freedom and unadjusted exploratory p value.
- `exploratory_tau_squared_dl`, `exploratory_tau_dl`: DerSimonian–Laird between-cohort variance and standard deviation on the stated scale.
- `effect_scale`, `se_method`, `assumption`, `interpretation`: scale, standard-error approximation and inferential limitations. The figure uses `logit_auc`.

Percentile-CI-derived standard errors and overlapping LODO training sets limit inference; these are exploratory annotations, not confirmatory heterogeneity tests.
