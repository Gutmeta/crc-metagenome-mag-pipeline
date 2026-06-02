# Workflow Overview

This repository preserves a representative set of numbered analysis scripts used in the CRC metagenome MAG workflow. Parameter-trial scripts and repeated cohort-specific scripts are collapsed to one example per step.

## 1. Read Processing, Assembly, and MAG Generation

- `01_` to `03_`: generate and run per-sample preprocessing, assembly/binning, and DiTASiC setup jobs.
- `04_` to `05_`: run GTDB-Tk classification and merge abundance tables with taxonomy.

## 2. Co-abundance Network Analysis

- `06_`: prepare group-specific FastSpar input tables using the representative core MAG prevalence filter.
- `07_`: run FastSpar correlation and bootstrap jobs.
- `08_` to `11_`: extract significant edges, summarize stable networks, plot correlation patterns, and run WGCNA.

## 3. CCDC1 / C-TCG Derivation

- `13_` to `18_`: extract high-quality MAGs, run CCDC1 DiTASiC workflows, and provide one representative ROC evaluation script.
- `19_` to `27_`: provide representative abundance, network, plotting, and WGCNA scripts for validation datasets.
- `28_` to `31_`: classify C1A/C1B clusters, provide one representative cluster reassignment/ANI workflow, and dereplicate genomes.
- `32_` to `38_`: extract C-TCG genomes, quantify them with DiTASiC, rank features, and extract selected representative FASTA files.

## 4. CC-TCG Construction and Validation

- `39_` to `40_`: run CC-TCG DiTASiC workflows, GTDB-Tk tree/classification, iTOL annotation generation, and abundance matrix construction.
- `41_` to `42_`: provide representative external validation DiTASiC and ROC evaluation examples.

## Notes for Reuse

- Replace `/path/to/...` placeholders with local paths or convert the relevant scripts to read from `config/example_config.yaml`.
- Treat files ending in `_example` as cohort templates; copy and adapt them for each dataset.
- Update Slurm settings before running on a different cluster.
- Keep raw data and generated outputs outside the Git repository.
