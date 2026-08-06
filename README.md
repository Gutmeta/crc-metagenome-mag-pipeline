# CRC Metagenome MAG Pipeline

Code and documentation for a colorectal cancer (CRC) metagenome and
metagenome-assembled genome (MAG) analysis workflow.

This repository accompanies a manuscript describing read processing, MAG
profiling, taxonomic annotation, co-abundance network analysis, C-TCG/CC-TCG
feature construction, external validation, and downstream machine-learning
evaluation. The numbered scripts and plotting interfaces document the
computational procedures used at each stage.

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

## Repository Contents

- Numbered scripts covering preprocessing, MAG generation and profiling,
  taxonomic annotation, co-abundance networks, C-TCG/CC-TCG construction, and
  validation analyses.
- Configuration and sample-metadata templates for adapting the workflow to a
  local compute environment.
- A stage-by-stage [workflow overview](docs/workflow_overview.md) that maps the
  numbered scripts to their analytical roles.
- Thirteen command-line plotting examples for networks, association summaries,
  model evaluation, annotated phylogenies, culture representation, ordination,
  and functional evidence.
- A Conda environment specification for the shared Python, R, and command-line
  dependencies.

## Data Availability and Repository Scope

This repository distributes source code, configuration templates, and a
synthetic plotting-data generator. Study sequencing data and derived
participant-level results are not bundled with the code.

- Study inputs should be obtained through the repositories, accession records,
  and data-access routes reported in the accompanying manuscript.
- GTDB, host-filtering, and other reference databases must be downloaded from
  their respective providers under the applicable terms.
- Large intermediate files, trained-model outputs, and generated figures should
  be stored outside the Git repository.
- Filesystem locations are expressed as `/path/to/...` placeholders and must be
  configured for the user's compute environment.

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
[`docs/plotting_examples.md`](docs/plotting_examples.md). To run a synthetic
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

Execution of the numbered workflow requires local configuration of input paths,
reference databases, software environments, and cluster resources. These
settings are kept explicit so that the analytical steps can be audited and
adapted to the target compute environment.

## Reproducibility

- `environment.yml` records the shared Python, R, and command-line dependencies.
- The numbered scripts correspond to the principal analytical stages described
  in the manuscript; dataset-specific variants use the same documented workflow
  interfaces.
- The plotting examples demonstrate figure generation with synthetic inputs;
  manuscript values must be supplied from the corresponding analysis outputs.
- Record the repository URL and exact release tag or commit identifier with each
  analysis so that the software version remains traceable.

The repository is released under the [MIT License](LICENSE).

## Citation

To cite this software, use the repository title, URL, and exact release tag or
commit identifier. For analyses reported in an associated manuscript or
preprint, cite the corresponding article in addition to the software version.
