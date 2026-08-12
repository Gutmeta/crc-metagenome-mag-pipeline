# Reference results

This directory contains public-safe snapshots of the current prediction, robustness, and final-figure outputs. It intentionally excludes `figure1_predictions.tsv`, `figure1_final_signature_predictions.tsv`, and other sample-level outputs carrying source aliases.

- `prediction/`: nested LODO and locked-panel cohort summaries, feature lists, selection curves, and the integrated result figure.
- `robustness/`: aggregate AUC comparisons, Top-N and named-species audits, covariate-adjusted disease rankings, and claim-support reports.
- `final_figure/`: CRC/CRA-merged overlap figure and its aggregate source tables.

The association and species-overlap reference outputs are kept beside their generating scripts in `key_species_overlap_v2/` and `species_association_overlap/`, because those scripts overwrite the same files during reproduction.

Prediction tables distinguish the strictly nested analysis from locked-panel estimates. Feature identities for locked panels were selected using the full disease dataset; their cohort AUCs are post-selection estimates.
