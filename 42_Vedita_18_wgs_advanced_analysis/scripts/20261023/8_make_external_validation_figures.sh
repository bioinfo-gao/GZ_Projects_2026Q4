#!/usr/bin/env bash
# P42 stage 8: publication-style external-validation audit figures.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/0_common.sh"

STAGE7="$OUT/07_external_sv_dld1_normalization_v1"
FIGDIR="$PROJECT/figures/20261023_fullcohort_repeat_resolve_v12"
[[ -s "$STAGE7/.complete" ]] || { echo "Stage 7 is incomplete: $STAGE7" >&2; exit 2; }
[[ ! -e "$FIGDIR" ]] || { echo "Refusing to overwrite figure directory: $FIGDIR" >&2; exit 2; }
"$PY" "$SCRIPT_DIR/8_make_external_validation_figures.py" \
  --stage7 "$STAGE7" --sample-table "$SAMPLES" --out-dir "$FIGDIR"
[[ -s "$FIGDIR/figure_manifest.json" && -s "$FIGDIR/P42_Fig1_external_caller_audit.png" \
  && -s "$FIGDIR/P42_Fig2_DLD1_GRIDSS_support_matrix.png" && -s "$FIGDIR/P42_Fig3_full_evidence_cluster_audit.png" ]] || {
  echo "Figure stage incomplete" >&2; exit 2;
}
