# CRC Metagenome MAG Pipeline

Example workflow scripts for a colorectal cancer (CRC) metagenome and MAG-based analysis pipeline.

This repository is intended as companion code for a manuscript. It provides a representative set of numbered workflow scripts used for read processing, MAG profiling, taxonomic annotation, co-abundance network analysis, C-TCG/CC-TCG feature construction, and downstream machine-learning evaluation.

## Repository Layout

```text
crc-metagenome-mag-pipeline/
├── README.md
├── LICENSE
├── environment.yml
├── config/
│   └── example_config.yaml
├── metadata/
│   └── example_sample_metadata.tsv
├── scripts/
│   ├── 01_... to 42_...
│   └── plotting/
│       └── portable plotting examples
└── docs/
    ├── workflow_overview.md
    └── plotting_examples.md
```

## What Is Included

- Representative workflow scripts from steps 1-42, renamed with two-digit prefixes.
- One default example for parameter trials, such as `06_prepare_fastspar_input.py` for the core MAG prevalence filter.
- One representative `*_example` script for repeated dataset-specific analyses.
- Example configuration and sample metadata templates.
- A compact workflow overview for readers and reviewers.
- Twelve privacy-safe plotting examples covering networks, heatmaps, associations,
  model summaries, annotated phylogenies, ordination, and functional evidence.

## What Is Not Included

- Raw sequencing files, SRA files, FASTQ files, BAM files, MAG FASTA files, or large intermediate results.
- Auto-generated per-sample job directories such as `1.jobs/`, `3.DiTASiC_jobs/`, `17.zkzzfx.jobs/`, and `19.DiTASiC_jobs/`.
- Redundant trial scripts for alternative thresholds or cohort-specific repetitions. These were collapsed to one representative example per workflow step.
- Machine-specific private paths. Original local paths were replaced with placeholders such as `/path/to/data2`, `/path/to/storage`, `/path/to/conda`, and `/path/to/crc-metagenome-mag-pipeline`.

## Dependencies

Create the base Conda environment with:

```bash
conda env create -f environment.yml
conda activate crc-mag-pipeline
```

Several workflow steps also require external command-line tools and reference databases, including SRA Toolkit, KneadData, Trimmomatic, MEGAHIT/metaSPAdes-style assemblers depending on the wrapper used, CheckM, dRep, GTDB-Tk, FastSpar, FastANI, DiTASiC, Kallisto, SeqKit, and GTDB/host reference databases.

## Usage

1. Edit `config/example_config.yaml` for your local data, database, Conda, and tool paths.
2. Prepare a sample metadata table following `metadata/example_sample_metadata.tsv`.
3. Run the numbered scripts in order for the analysis section you need. For numbered `*_example` scripts, edit the dataset name and placeholder paths before use.
4. For Slurm-based scripts, update `#SBATCH` resources, partitions, and log directories before submission.

Scripts under `scripts/plotting/` use command-line arguments and do not require
editing source paths. Their complete input schemas are documented in
[`docs/plotting_examples.md`](docs/plotting_examples.md). To run a privacy-safe
smoke test entirely outside the repository:

```bash
python scripts/plotting/generate_demo_inputs.py \
  --output-dir /tmp/crc-mag-plotting-demo/inputs
mkdir -p /tmp/crc-mag-plotting-demo/figures

python scripts/plotting/plot_association_heatmap_example.py \
  --input /tmp/crc-mag-plotting-demo/inputs/association_heatmap.tsv \
  --output-prefix /tmp/crc-mag-plotting-demo/figures/association_heatmap
```

The example writes editable PDF/SVG figures and a PNG preview. All generated
records use deterministic fictitious identifiers such as `demo_001`.

Numbered workflow example:

```bash
bash scripts/02_pipe.sh <sample_list.txt> <input_fastq_root> <output_dir>
python scripts/06_prepare_fastspar_input.py
bash scripts/07_run_fastspar_1.sh
python scripts/08_extract_network_edges.py
```

The scripts are presented as a transparent workflow record for publication. They are not guaranteed to be one-command portable without adapting paths, databases, cluster settings, and dataset metadata.

## Reproducibility and Data Scope

- `environment.yml` records the shared Python, R, and command-line dependencies.
- The numbered scripts preserve representative manuscript analysis logic; repeated cohort-specific and parameter-trial scripts are intentionally collapsed to one documented example.
- The plotting examples reproduce figure-generation patterns with synthetic inputs, not the manuscript's numerical results or participant-level analyses.
- Raw data, private metadata, model predictions, generated figures, and machine-specific paths are not distributed in this repository.
- For peer review or publication, report the repository URL together with a release tag or commit identifier so the reviewed code remains unambiguous.

The repository is released under the [MIT License](LICENSE).

## Citation

For peer review, cite the accompanying manuscript and include the repository URL
and reviewed commit or release. Replace this note with the manuscript DOI or
preprint citation when a persistent identifier becomes available.
