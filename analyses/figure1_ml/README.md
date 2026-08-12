# Figure 1 ML reproducibility package

This package contains portable analysis code, minimized pseudonymized analysis inputs, English methods, provenance records, reference results, and integrity checks for the species-level microbiome analyses of colorectal adenoma (CRA), colorectal cancer (CRC), inflammatory bowel disease (IBD), and irritable bowel syndrome (IBS).

## Released inputs and privacy boundary

This directory is the public reproducibility package for the Figure 1 ML
analyses. Machine-specific paths, direct source sample identifiers, coordinates,
dates, free-text sample descriptions, medication fields, and unused clinical
variables were removed from the analysis inputs. Retained sample-level fields
are limited to disease status, age, sex, and taxonomic relative abundance.

Pseudonymization is not irreversible anonymization. Age, sex, disease status,
and high-dimensional abundance profiles can permit linkage when combined,
particularly when source records are public. The release therefore documents
the residual risk, retained fields, provenance, and applicable data terms in
`DEIDENTIFICATION.md`, `DATA_SOURCES.md`, and `LICENSE_NOTICE.md`.

## Included data

The offline analysis inputs comprise:

- 20 minimized metadata files: 6,348 rows; and
- 20 MetaPhlAn 4 species files: 6,264 profiled samples and 934,240 long-format rows.

Every included input is listed with its byte count and SHA-256 digest in `DATA_MANIFEST.tsv`. Original Metalog exports are not included. No source-to-pseudonym crosswalk is included or written by the release builder.

## Directory layout

- `CRA/`, `CRC/`, `IBD/`, and `IBS/`: minimized inputs arranged by study; Metalog acquisition scripts are retained for provenance and optional reacquisition.
- `key_species_overlap_v2/`: cohort screening, association testing, random-forest stability ranking, meta-analysis, and reproduced outputs.
- `requested_topjournal_figures/`: nested leave-one-dataset-out modelling and figure code.
- `species_association_overlap/`: Top-N association-set construction, cross-disease permutation tests, standardized IBS association inputs, and reproduced outputs.
- `topjournal_robustness_v1_20260530/`: AUC comparisons, age/sex-adjusted sensitivity analyses, and final merged-figure code.
- `reference_results/`: public-safe prediction, robustness, and final-figure tables; sample-level prediction files with original identifiers are intentionally excluded.
- `provenance/`: sanitized source manifests, retrieval/validation timestamps, cohort inclusion/exclusion records, and source inventory.
- `tools/`: input-release builder and public-release validator.

Rendered PNG and PDF files are intentionally not versioned. They are generated
by the commands below from the included inputs and aggregate tables.

## Quick validation and reproduction

Use Python 3.10 or later. Install the frozen environment, validate the release, and run the offline core workflow:

```bash
python3 -m venv .venv
.venv/bin/pip install -r requirements.txt
.venv/bin/python tools/validate_public_release.py
bash run_pipeline.sh core .venv/bin/python
```

The core workflow runs the primary MetaPhlAn 4 association analysis, species-overlap analysis, and overlap figures. The full workflow additionally reruns nested LODO prediction, robustness analyses, and the final merged figure:

```bash
bash run_pipeline.sh all .venv/bin/python
```

The full LODO and robustness workflow is computationally more expensive. It
recreates aggregate analysis outputs and untracked figures within the package
directories but does not modify the packaged inputs.

## Reacquiring source exports

The release inputs are already included, so downloading is not required to reproduce the package. The acquisition scripts use Metalog `latest` URLs and are retained only for source tracing or rebuilding from current upstream exports. A future download can differ from the version inspected on 24 May 2026 and can overwrite release-formatted files; run acquisition in a separate working copy.

## Validation

The minimized inputs reproduced the primary association/ranking outputs and all 11 species-association overlap outputs byte-for-byte relative to the working analysis. The final validation report is `VALIDATION_REPORT.txt`. Random seeds and thresholds remain fixed in the scripts.

## Statistical reporting note

`run_lodo_bestk_auc.py` writes both strictly nested LODO estimates and post-selection locked-panel estimates. Only the nested analysis keeps the outer held-out cohort outside feature selection and supports unbiased generalization claims. Locked-panel cohort AUCs are post-selection performance estimates.

The complete manuscript-ready English methods are in `METHODS.md`; variable definitions are in `DATA_DICTIONARY.md`; source timing is in `DATA_SOURCES.md`.

## Licensing

Original analysis code is covered by the repository's MIT License. Metalog's
database and database contents are distributed under the ODbL 1.0 and DbCL 1.0,
respectively. External aggregate association tables retain their source terms
and citation requirements. See `LICENSE_NOTICE.md` and `LICENSES/README.md`.
