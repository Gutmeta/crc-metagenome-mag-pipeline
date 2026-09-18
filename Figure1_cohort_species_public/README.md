# Figure 1: frozen cohort-level species abundance inputs

This companion dataset contains the **20 independent cohort species tables** underlying the final Figure 1 analysis: 10 colorectal cohorts, 9 IBD cohorts and one IBS cohort (Mars_2020). Some colorectal studies contribute both CRA and CRC comparisons; they are not duplicated as independent source tables.

## Contents

- `profiles/`: 20 compressed long-format MetaPhlAn 4 species relative-abundance tables.
- `sample_groups/`: 20 minimal grouping tables containing a release sample key and applicable case/control labels only.
- `COHORT_TABLES.tsv` / `COHORT_TABLES.xlsx`: cohort names, established project accessions, table paths, original Metalog URLs, sizes and hashes.
- `ANALYSIS_GROUP_COUNTS.tsv`: case/control counts and file mapping for every disease comparison, including combined CRC/CRA.
- `VALIDATION_SUMMARY.json`: traceability and content-check summary.
- `MANIFEST.tsv`: frozen file sizes and SHA-256 digests.
- `scripts/validate.py`: checks file integrity, retained fields and analysis counts.
- `scripts/export_matrix.py`: optional conversion to species-by-sample matrices.

The package contains **6,264 profiled samples and 934,240 abundance rows**. Samples are not necessarily distinct participants: repeated visits or longitudinal specimens can occur. Disease comparisons reuse some controls and must not be summed to infer independent sample counts.

## What is original, and what has been changed?

The original downloaded Metalog species exports were separately locked byte-for-byte for local provenance, together with the exact de-identified upstream analysis inputs. The original export source URLs and SHA-256 hashes appear in `COHORT_TABLES.tsv`.

This public companion uses the exact pseudonymized, analysis-eligible species profiles from the existing upstream release inputs. Its profile files are copied without rewriting. Every retained `(sample, species, abundance)` record was compared with the original export after mapping sample keys in memory; all retained taxon names and abundance strings matched exactly. No additional Top-300 selection, prevalence threshold, numerical rounding or abundance renormalization has been applied here.

Original downloaded exports collectively contain 6,303 profiled samples. The public analysis-input subset contains 6,264 because the upstream release excludes samples outside the eligible disease comparisons. Therefore the public files must be described as **pseudonymized analysis-eligible species profiles**, not untouched full original exports. Source-to-release identifier mappings are not included or written.

The table contains microbial species **relative abundances**, not raw sequencing reads or absolute organism counts. Values are retained in their source-export scale. Use the source values as supplied; do not silently rescale them.

## Table schema

Each `profiles/*.species.tsv.gz` has three columns:

| Column | Meaning |
|---|---|
| `sample_alias` | Release identifier matching `sample_groups/*.groups.tsv` column `sample_id` |
| `clade_name` | Species identifier starting with `s__`, including source SGB names |
| `rel_abund` | Original exported relative-abundance value |

Grouping files contain `sample_id` and applicable `CRA`, `CRC`, `CRC_CRA`, `IBD` or `IBS` columns. Values are `case`, `control`, or `not_in_comparison`. For example, a colorectal carcinoma sample can be a case for CRC and CRC_CRA while not belonging to the CRA comparison. The broader diagnosis strings, age, sex, dates, geography and free-text fields are not included.

## Traceability to the final figure

The identified upstream package was checked against the final figure inputs:

- All 24 AUC estimates and their CI endpoints, held-out case/control counts and selected K matched.
- All 681 selected-set membership records matched.
- All 3,635 detectable-background membership records matched.
- All 20 profile files matched the upstream release's recorded SHA-256 hashes.
- Original-export abundance values matched every retained pseudonymized row.

This establishes table-level provenance. It does not constitute a new rerun of all nested model training. The figure-level code package remains separate; adding these profiles does not automatically turn its plotting code into a complete raw-profile model-training pipeline.

The Mars_2020 project accession was not recorded in the supplied Figure 1 accession map; this dataset retains the public study code and original Metalog source URL instead of guessing an accession.

## Source version and public-data scope

The source URLs contain `latest` and may change. The frozen SHA-256 digests, not a future download, identify the files used here. Recorded acquisition/validation times come from the original source manifests; a reuse/validation timestamp does not establish the exact original transfer time. A missing timestamp means the supplied manifest did not record one.

The original upstream package describes its sample-level release as a pseudonymized public-release candidate, subject to the source-data redistribution terms. This preparation verifies provenance and removes unused identifying fields; it does not establish new permissions to redistribute another provider's data. Pseudonymization is not a guarantee of irreversible anonymization. Do not add the local raw-export snapshot, internal path registry or identifier crosswalk to a public repository.

## Use

Install pandas and NumPy, then validate:

```bash
python scripts/validate.py
```

To produce a species-by-sample matrix in a separate location:

```bash
python scripts/export_matrix.py --cohort PRJEB7774 --output ../PRJEB7774_species_matrix.tsv.gz
```

For a matrix, missing sample-species pairs are filled with zero and repeated pairs are summed, matching the upstream loading convention. The frozen long-format files remain unchanged. Keep this whole folder together when sharing it as a companion data deposit.
