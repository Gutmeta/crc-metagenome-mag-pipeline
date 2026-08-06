# Portable Plotting Examples

The scripts in `scripts/plotting/` preserve representative plotting patterns from the CRC, IBD, pan-disease, and TCG-FDC analyses while removing project-specific paths, private identifiers, fixed internal batch names, and patient-level data.

Each visual type is represented once across the source projects. The examples consume small, documented tables and write editable PDF/SVG figures plus optional PNG previews. They do not train predictive models or require access to the manuscript data.

## Quick Start

Create the Conda environment described in the repository root, then generate deterministic synthetic inputs outside the repository:

```bash
conda env create -f environment.yml
conda activate crc-mag-pipeline

python scripts/plotting/generate_demo_inputs.py \
  --output-dir /tmp/crc-mag-plotting-demo/inputs
mkdir -p /tmp/crc-mag-plotting-demo/figures
```

The generated records use identifiers such as `demo_001`; they are entirely synthetic and contain no study participant information.

All plotting scripts accept `--formats pdf,svg,png`. Use a subset such as `--formats pdf` when only one output is needed.

## Included Figure Types

| Example | Representative lineage | Input |
|---|---|---|
| `plot_consensus_network_example.py` | `pipeline_CRC` consensus networks | Node and edge tables |
| `plot_association_heatmap_example.py` | `pipeline_IBD` clinical heatmaps | Long-form effect/FDR table |
| `plot_association_scatter_example.py` | `pipeline_IBD` within-cohort associations | Grouped x/y observations |
| `plot_performance_matrix_example.py` | `pipeline_Pan` model comparisons | Paired AUROC table |
| `plot_annotated_phylogeny_example.R` | `pipeline_Pan` annotated trees | Newick tree and tip annotations |
| `plot_group_distributions_example.R` | `pipeline_Pan` case-control shifts | Cohort/metric/group observations |
| `plot_volcano_example.R` | `pipeline_Pan` cross-cohort summaries | Precomputed effects and FDR values |
| `plot_auc_forest_example.py` | `TCG-FDC` AUROC summaries | Precomputed AUROC confidence intervals |
| `plot_roc_curves_example.py` | `TCG-FDC` ROC panels | Out-of-fold or external prediction scores |
| `plot_upset_example.py` | `TCG-FDC` structure/function overlap | Long-form set membership |
| `plot_pcoa_example.py` | `TCG-FDC`/`pipeline_IBD` ordination | Square distance matrix and metadata |
| `plot_evidence_bubbles_example.py` | `TCG-FDC` functional summaries | Long-form effect/FDR table |

## Commands and Input Schemas

### Consensus network

The node table requires `node_id`, `x`, `y`, `group`, and positive `size` columns. The edge table requires `source`, `target`, `weight`, `sign` (`positive` or `negative`), and `support`. Every edge endpoint must occur in the node table.

```bash
python scripts/plotting/plot_consensus_network_example.py \
  --nodes /tmp/crc-mag-plotting-demo/inputs/network_nodes.tsv \
  --edges /tmp/crc-mag-plotting-demo/inputs/network_edges.tsv \
  --output-prefix /tmp/crc-mag-plotting-demo/figures/consensus_network
```

### Association heatmap

The long-form table requires `row`, `column`, `effect`, and `q_value`. Row/column pairs must be unique, and FDR values must lie in `[0, 1]`.

```bash
python scripts/plotting/plot_association_heatmap_example.py \
  --input /tmp/crc-mag-plotting-demo/inputs/association_heatmap.tsv \
  --output-prefix /tmp/crc-mag-plotting-demo/figures/association_heatmap
```

### Association scatter plot

The table requires numeric `x` and `y` values plus a `group`. Each group needs at least four observations and variation in `x`. Confidence bands are ordinary linear-model intervals for visualization; repeated-measures or cohort-specific statistics should be computed upstream.

```bash
python scripts/plotting/plot_association_scatter_example.py \
  --input /tmp/crc-mag-plotting-demo/inputs/association_scatter.tsv \
  --x-label "Guild score" --y-label "Clinical measure" \
  --output-prefix /tmp/crc-mag-plotting-demo/figures/association_scatter
```

### Performance comparison matrix

The table requires unique `dataset` values and `model_a_auc`/`model_b_auc` values in `[0, 1]`. The displayed delta is calculated as model B minus model A.

```bash
python scripts/plotting/plot_performance_matrix_example.py \
  --input /tmp/crc-mag-plotting-demo/inputs/performance_matrix.tsv \
  --output-prefix /tmp/crc-mag-plotting-demo/figures/performance_matrix
```

### Annotated phylogeny

The annotation table requires `tip_id`, `group`, `status`, and a finite, non-negative `score`. Every Newick tip must have one annotation row.

```bash
Rscript scripts/plotting/plot_annotated_phylogeny_example.R \
  --tree /tmp/crc-mag-plotting-demo/inputs/demo_tree.nwk \
  --annotation /tmp/crc-mag-plotting-demo/inputs/tree_annotation.tsv \
  --output-prefix /tmp/crc-mag-plotting-demo/figures/annotated_phylogeny
```

### Group distributions

The table requires `cohort`, `metric`, `group`, and numeric `value`. The portable example expects exactly two group labels and calculates a Wilcoxon test within each cohort/metric panel, followed by Benjamini-Hochberg correction across panels.

```bash
Rscript scripts/plotting/plot_group_distributions_example.R \
  --input /tmp/crc-mag-plotting-demo/inputs/group_distributions.tsv \
  --output-prefix /tmp/crc-mag-plotting-demo/figures/group_distributions
```

### Volcano plot

The table requires unique `feature` values, precomputed `effect`, and `q_value` in `(0, 1]`. An optional `label` column controls preferred annotations; otherwise the most significant features are labeled.

```bash
Rscript scripts/plotting/plot_volcano_example.R \
  --input /tmp/crc-mag-plotting-demo/inputs/volcano.tsv \
  --effect-threshold 0.5 --q-threshold 0.05 \
  --output-prefix /tmp/crc-mag-plotting-demo/figures/volcano
```

### AUROC forest plot

The table requires `comparison`, `estimate`, `ci_low`, `ci_high`, `n_positive`, and `n_negative`. Each confidence interval must satisfy `0 <= ci_low <= estimate <= ci_high <= 1`.

```bash
python scripts/plotting/plot_auc_forest_example.py \
  --input /tmp/crc-mag-plotting-demo/inputs/auc_forest.tsv \
  --output-prefix /tmp/crc-mag-plotting-demo/figures/auc_forest
```

### ROC curves

The prediction table requires `comparison`, binary `label`, and a probability-like `score` in `[0, 1]`. Each comparison must contain both classes. The script calculates AUROC and a class-stratified bootstrap confidence interval.

```bash
python scripts/plotting/plot_roc_curves_example.py \
  --input /tmp/crc-mag-plotting-demo/inputs/roc_predictions.tsv \
  --bootstrap 1000 --seed 42 \
  --output-prefix /tmp/crc-mag-plotting-demo/figures/roc_curves
```

### UpSet-style intersections

The long-form table requires `item_id` and `set_name`, with one row per membership. Duplicate membership rows are ignored.

```bash
python scripts/plotting/plot_upset_example.py \
  --input /tmp/crc-mag-plotting-demo/inputs/set_membership.tsv \
  --output-prefix /tmp/crc-mag-plotting-demo/figures/upset
```

### PCoA

The distance table must be square: its first column contains sample identifiers and the remaining column names contain the same identifiers. Values must be finite, non-negative, symmetric, and zero on the diagonal. Metadata requires unique `sample_id` and `group` values.

```bash
python scripts/plotting/plot_pcoa_example.py \
  --distance-matrix /tmp/crc-mag-plotting-demo/inputs/distance_matrix.tsv \
  --metadata /tmp/crc-mag-plotting-demo/inputs/ordination_metadata.tsv \
  --output-prefix /tmp/crc-mag-plotting-demo/figures/pcoa
```

### Functional evidence bubbles

The table requires unique `feature`/`database`/`comparison` combinations, numeric `effect`, and `q_value` in `(0, 1]`. Color encodes effect direction and size encodes `-log10(q)`.

```bash
python scripts/plotting/plot_evidence_bubbles_example.py \
  --input /tmp/crc-mag-plotting-demo/inputs/evidence_bubbles.tsv \
  --output-prefix /tmp/crc-mag-plotting-demo/figures/evidence_bubbles
```

## Privacy and Reuse Notes

- Supply only publication-safe, de-identified tables. The plotting scripts do not need names, email addresses, clinical record numbers, cluster paths, or raw sequencing identifiers.
- Keep generated tables and figures outside the repository. Existing `.gitignore` rules exclude common tabular and image outputs.
- The examples intentionally separate visualization from cohort-specific modeling. Reproduce the statistical design appropriate to a study before passing estimates or predictions into these scripts.
- PDF and SVG outputs keep text editable when the installed graphics backend supports it; fonts fall back through Arial, Liberation Sans, and DejaVu Sans.
