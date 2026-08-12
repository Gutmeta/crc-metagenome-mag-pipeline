# Data dictionary

## Minimized metadata: `*/<study>/metadata_all_wide.tsv`

| Field | Type | Meaning | Release handling |
|---|---|---|---|
| `study_code` | string | Public study identifier used by Metalog and the analysis | Retained for cohort-level validation and provenance |
| `sample_alias` | string | Release pseudonym in the form `SMP_` plus 20 hexadecimal characters | Replaces the source sample identifier consistently within the metadata and profile files; no crosswalk is released |
| `subject_disease_status` | string | Source disease-status category used to define cases and controls | Retained without relabelling so the published case/control rules can be audited |
| `age_years` | numeric or missing | Age in years as supplied by the source | Retained only for the age/sex-adjusted sensitivity analysis |
| `sex` | categorical or missing | Sex category as supplied by the source | Retained only for the age/sex-adjusted sensitivity analysis |

Metadata files are tab-separated UTF-8 text with one row per retained sample within a study. Missing age or sex is represented by an empty field. Exact disease-label mappings and eligibility thresholds are defined in `METHODS.md` and in the analysis scripts.

## MetaPhlAn 4 species profiles: `metaphlan4_species_*_latest.tsv.gz`

| Field | Type | Meaning |
|---|---|---|
| `sample_alias` | string | Pseudonym matching the study metadata |
| `clade_name` | string | MetaPhlAn species-level clade; only values beginning with `s__` are retained |
| `rel_abund` | numeric | Source-provided relative abundance; values are retained without rescaling |

Files are gzip-compressed long-format TSVs. Absent sample-species pairs are interpreted as zero when the analysis pivots the file to a matrix.

## Manifests and result tables

- `DATA_MANIFEST.tsv` is the authoritative inventory of released sample-level inputs.
- `provenance/disease_arm_inventory.tsv` records candidate arms, inclusion decisions, and aggregate case/control counts; its paths are package-relative.
- `key_species_overlap_v2/` and `species_association_overlap/` contain reproduced aggregate result tables.
- `reference_results/` contains prediction and robustness tables that do not expose original sample aliases. Original sample-level prediction tables are deliberately omitted.
