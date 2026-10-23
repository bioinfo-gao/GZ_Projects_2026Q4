#!/usr/bin/env bash
# P42 stage 10: all client-delivery figures from completed audited stages.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/0_common.sh"
FIGDIR="$PROJECT/figures/20261023_fullcohort_repeat_resolve_delivery_v7"

[[ -s "$OUT/09_fixed_candidate_regenotyping_v1/.complete" ]] || { echo "Stage 9 is incomplete" >&2; exit 2; }
[[ ! -e "$FIGDIR" ]] || { echo "Refusing to overwrite figure directory: $FIGDIR" >&2; exit 2; }
"$PY" "$SCRIPT_DIR/10_make_delivery_figure_suite.py" --project "$PROJECT" --out-dir "$FIGDIR"
[[ -s "$FIGDIR/figure_manifest.json" ]] || { echo "Missing figure manifest" >&2; exit 2; }
shopt -s nullglob
pngs=("$FIGDIR"/P42_Fig*.png)
[[ "${#pngs[@]}" -eq 12 ]] || { echo "Expected 12 PNG figures" >&2; exit 2; }
