# Provenance files

- `CRA_download_manifest.tsv`: Metalog URLs, original acquisition status, and recorded byte counts for the retained CRA metadata and MetaPhlAn 4 species files. This manifest drives `CRA/download_metalog_cra_from_manifest.py`.
- `CRC_metalog_download_manifest.tsv`, `IBD_metalog_download_manifest.tsv`, and `IBS_metalog_download_manifest.tsv`: metadata and MetaPhlAn 4 species acquisition or validation records, including the local timestamp recorded by the original scripts on 24 May 2026. A status of `existing_local_file` or `already_exists` indicates validation or reuse rather than a new network transfer at the recorded time.
- `metalog_gut_disease_inventory.tsv`: MetaPhlAn 4 species availability summary across candidate disease studies.
- `disease_arm_inventory.tsv`: auditable inclusion table for every candidate disease arm, including matched case/control counts and exclusion reasons.
- `figure1_exclusions.tsv`: loading or eligibility exclusions produced by the prediction workflow.

Machine-specific path prefixes were removed from these public-package copies. Relative paths are resolved from the package root.

The files in this directory describe the original source audit and are not the released sample-level inputs. The public sample-level inputs are the minimized files under `CRA/`, `CRC/`, `IBD/`, and `IBS/`, inventoried in `DATA_MANIFEST.tsv`. Source byte counts in the acquisition manifests therefore differ from the minimized-release byte counts.
