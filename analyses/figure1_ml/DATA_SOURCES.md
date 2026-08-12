# Data sources and timestamps

## Metalog inputs

Human metadata and precomputed taxonomic profiles were obtained from Metalog (`https://metalog.embl.de/`). Source endpoints for every inspected package are recorded in `provenance/`.

Metalog identifies its database as ODbL 1.0 and its individual database
contents as DbCL 1.0. Reusers must retain Metalog attribution and comply with
the applicable share-alike and content-license terms. License scope and links
are recorded in `LICENSE_NOTICE.md` and `LICENSES/README.md`.

The timestamps below are local Asia/Shanghai time (`UTC+08:00`) recorded by the original scripts or filesystem. A manifest status of `existing_local_file` or `already_exists` means that the file was validated or reused at that time; it does **not** prove that the network transfer occurred at that exact timestamp.

| Source group | Recorded interval | Interpretation |
|---|---|---|
| IBD | 2026-05-24 13:02:52–13:07:57 | Retained metadata and MetaPhlAn 4 retrieval/validation records |
| IBS | 2026-05-24 14:11:37–14:11:43 | Retained metadata and MetaPhlAn 4 retrieval/validation records |
| CRC | 2026-05-24 14:41:56–14:43:05 | Retained metadata and MetaPhlAn 4 retrieval/validation records |
| CRA `Lee_2023_CRC` metadata | 2026-05-24 11:55:11 | Original local file modification time; CRA manifest has no timestamp column |
| CRA `Lee_2023_CRC` MetaPhlAn species profile | 2026-05-24 11:55:22 | Original local file modification time; CRA manifest has no timestamp column |

Therefore, the defensible summary is: the principal source acquisition/validation session occurred on **24 May 2026**, but the exact download time cannot be established for every file from the retained evidence. The manifests are the authoritative audit records.

Metalog profile URLs contain `latest`. Reacquisition at a later date is not guaranteed to return byte-identical input. The release package therefore contains the minimized analysis inputs and their SHA-256 hashes.

## External IBS association inputs

Two standardized IBS association tables used to augment the single primary Metalog IBS arm are included under `species_association_overlap/external_ibs_associations/`:

- `Goll_2020_FMT_IBS_baseline_vs_donor.standardized.tsv`; and
- `PRJEB34103_IBS_vs_HC_deseq2_species_s3.standardized.tsv`.

Their comparison direction, sample counts, and source notes are recorded in `external_ibs_association_source_summary.tsv`. These are aggregate species-association tables rather than participant-level metadata.

## Release derivation

The minimized package was assembled on 12 August 2026. `DATA_MANIFEST.tsv` records every released input and digest. `tools/build_deidentified_inputs.py` documents the transformation; original exports and identifier crosswalks are excluded.
