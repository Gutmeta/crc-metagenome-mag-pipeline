# Validation and reproduction scope

The repository provides figure-level inputs, IBS species-level selection evidence, numerical source data, supporting calculations, figure-reproduction code, public documentation and a ZIP archive of processed cohort abundance profiles with minimal disease-comparison labels. `MANIFEST.tsv` records distributed file sizes and SHA-256 hashes.

Input validation checks the AUC tables, selected-set sizes, background sizes and selected-species membership within each background. The IBS evidence table identifies all 465 background species and distinguishes the 81 random-forest candidates used in Figure 1 from their individual-species differential-test results. The companion archive contains its own file manifest, cohort/sample counts and field definitions.

Figure-level validation covers source-data consistency, workbook/TSV agreement, the 100,000 overlap draws, exploratory I² calculations and numerical statements in the figure legend. It also checks figure dimensions, Arial Regular typography, black text and plotting layout. The documented reproduction command regenerates figures, source tables, simulation counts, supporting calculations and manuscript descriptions from the released figure-level inputs. Comparisons use numerical results and image content because PDF, SVG and Office metadata can vary between runs.

These checks establish file integrity and consistency of the released figure-level calculations. They do not verify that every selected species is individually significant, demonstrate independent-participant sampling, or reproduce the initial nested model training and random-forest selection. The sample-level abundance archive and species-level IBS evidence extend the available data but do not supply the complete model-training pipeline. IBS remains a single-cohort candidate set.

See [README.md](../README.md) for the reproduction and validation commands, [Methods](Figure1_methods.md) for selection and statistical definitions, and [Data availability](Data_availability.md) for the released data scope.
