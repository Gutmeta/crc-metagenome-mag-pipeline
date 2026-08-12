# De-identification and data-minimization record

## Purpose and scope

This release transformation minimizes the Metalog-derived inputs to the variables required to rerun the documented analyses. It is a pseudonymization and data-minimization step, not a certification of anonymous data under any particular law or institutional policy.

## Retained information

- public study code;
- pseudonymized sample key;
- disease-status label;
- age in years and sex, when available; and
- species-level taxon names and abundance values.

Age and sex are retained because they are required for the prespecified covariate-adjusted sensitivity analysis. They are quasi-identifiers and should be removed if that sensitivity analysis is not part of a future release.

## Removed information

The release metadata exclude source sample aliases, subject identifiers, geographic coordinates, collection dates, raw sample titles and descriptions, medication fields, smoking/BMI/location fields, and all other unused clinical or free-text variables. MetaPhlAn rows above the species rank and samples not used in an included disease arm are also excluded.

## Pseudonym construction

For this release, a sample key is `SMP_` followed by the first 20 hexadecimal characters of a SHA-256 digest over a release namespace, study key, and source sample alias. The same key is used in the metadata and both abundance modalities within a study. The builder keeps the source-to-release mapping only in memory and never writes a crosswalk.

Because the construction is deterministic and source aliases may be public, it should not be treated as cryptographically irreversible anonymization. High-dimensional abundance profiles may also permit linkage to an upstream public record. The pseudonyms mainly prevent accidental disclosure and direct identifier propagation through analysis outputs.

## Rebuilding the release

`tools/build_deidentified_inputs.py` performs the transformation from a private/source working directory into a separate output root. It must never be run with a source root that points at the release package itself, and its output should be revalidated before publication.

## Release controls

This release retains exact age and sex because they are required for the
documented covariate-adjusted sensitivity analysis. The approved public scope
also includes disease status and species-level abundance profiles. Residual
linkage risk is disclosed above and should be considered when reusing or
combining the tables with other records.

Every release revision must pass `tools/validate_public_release.py` and the
checks recorded in `VALIDATION_REPORT.txt`. Source-to-pseudonym crosswalks,
source aliases, and sample-level prediction tables containing original
identifiers must never be added to the public package.
