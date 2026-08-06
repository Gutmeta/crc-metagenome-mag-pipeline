# Workflow Overview

The numbered scripts map the CRC metagenome and MAG analysis to four analytical
stages, from read processing through CC-TCG external validation. Dataset-specific
inputs, reference databases, and compute resources are supplied through the
documented configuration fields and script parameters.

## 1. Read Processing, Assembly, and MAG Generation

- `01_` to `03_`: generate and run per-sample preprocessing, assembly/binning, and DiTASiC setup jobs.
- `04_` to `05_`: run GTDB-Tk classification and merge abundance tables with taxonomy.

## 2. Co-abundance Network Analysis

- `06_`: prepare group-specific FastSpar input tables using the configured core-MAG prevalence filter.
- `07_`: run FastSpar correlation and bootstrap jobs.
- `08_` to `11_`: extract significant edges, summarize stable networks, plot correlation patterns, and run WGCNA.

## 3. CCDC1 / C-TCG Derivation

- `13_` to `18_`: extract high-quality MAGs, run CCDC1 DiTASiC workflows, and evaluate discrimination with ROC analysis.
- `19_` to `27_`: quantify validation cohorts and perform abundance, network, plotting, and WGCNA analyses.
- `28_` to `31_`: classify C1A/C1B clusters, evaluate cluster correspondence with ANI, and dereplicate genomes.
- `32_` to `38_`: extract C-TCG genomes, quantify them with DiTASiC, rank features, and export selected feature genomes.

## 4. CC-TCG Construction and Validation

- `39_` to `40_`: run CC-TCG DiTASiC workflows, GTDB-Tk tree/classification, iTOL annotation generation, and abundance matrix construction.
- `41_` to `42_`: quantify external validation cohorts with DiTASiC and evaluate ROC performance.

## Execution Notes

- Configure each `/path/to/...` placeholder with the corresponding local path,
  either directly or through `config/example_config.yaml` where supported.
- Files ending in `_example` expose reusable interfaces for analyses that require dataset-specific inputs or cohort definitions.
- Update Slurm settings before running on a different cluster.
- Keep raw data and generated outputs outside the Git repository.
