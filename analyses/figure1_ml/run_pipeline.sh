#!/usr/bin/env bash
set -eu

PACKAGE_ROOT=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
MODE=${1:-core}
PYTHON_BIN=${2:-python3}

case "$MODE" in
  core|all) ;;
  *)
    echo "Usage: bash run_pipeline.sh [core|all] [python-executable]" >&2
    exit 2
    ;;
esac

cd "$PACKAGE_ROOT"
export PYTHONDONTWRITEBYTECODE=1
MPLCONFIGDIR=$(mktemp -d "${TMPDIR:-/tmp}/figure1_ml_public_mpl.XXXXXX")
export MPLCONFIGDIR
trap 'rm -rf "$MPLCONFIGDIR"' EXIT
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-1}
export OPENBLAS_NUM_THREADS=${OPENBLAS_NUM_THREADS:-1}
export MKL_NUM_THREADS=${MKL_NUM_THREADS:-1}

"$PYTHON_BIN" tools/validate_public_release.py
"$PYTHON_BIN" key_species_overlap_v2/run_key_species_overlap_v2.py
"$PYTHON_BIN" species_association_overlap/run_species_association_overlap.py
"$PYTHON_BIN" species_association_overlap/make_publication_figures.py

if [ "$MODE" = "all" ]; then
  "$PYTHON_BIN" requested_topjournal_figures/run_lodo_bestk_auc.py
  "$PYTHON_BIN" requested_topjournal_figures/plot_requested_topjournal_figures.py
  "$PYTHON_BIN" requested_topjournal_figures/plot_feature_panel_auc_overlap_story.py
  "$PYTHON_BIN" topjournal_robustness_v1_20260530/scripts/run_topjournal_robustness_v1.py
  "$PYTHON_BIN" topjournal_robustness_v1_20260530/crc_cra_merged_overlap/scripts/make_crc_cra_merged_overlap_figure.py
fi

echo "Completed Figure1_ML $MODE workflow."
