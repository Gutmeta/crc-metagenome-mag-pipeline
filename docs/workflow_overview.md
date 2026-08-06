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

The culture-status helper consumes a query map, raw `skani search` results, and
GTDB metadata. By default it treats hits with ANI >=95% and maximum
bidirectional alignment fraction >=30% as species-level matches, assigns
`Cultured` when any qualifying matched species contains an isolate genome, and
writes both a tip-level evidence table and an iTOL color strip. Query genomes
remain MAGs; the status describes culture representation of the matched species.
The query map requires `tip_id` and `query_fasta`; the skani table uses its
standard `Ref_file`, `Query_file`, `ANI`, `Align_fraction_ref`, and
`Align_fraction_query` columns, with alignment fractions expressed as
percentages. GTDB metadata requires `accession`, `gtdb_taxonomy`, and
`ncbi_genome_category`.

```bash
python scripts/40_6_make_itol_cultured.py \
  --query-map query_genomes.tsv \
  --skani-results gtdb_representative_hits.tsv \
  --metadata bac120_metadata.tsv.gz \
  --status-output culture_status_evidence.tsv \
  --itol-output itol_culture_status.txt
```

## Execution Notes

- Configure each `/path/to/...` placeholder with the corresponding local path,
  either directly or through `config/example_config.yaml` where supported.
- Files ending in `_example` expose reusable interfaces for analyses that require dataset-specific inputs or cohort definitions.
- Update Slurm settings before running on a different cluster.
- Keep raw data and generated outputs outside the Git repository.
