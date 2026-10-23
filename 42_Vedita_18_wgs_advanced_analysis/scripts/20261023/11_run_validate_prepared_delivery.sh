#!/usr/bin/env bash
# P42 stage 11: fail-closed audit of the prepared but unsent delivery package.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/0_common.sh"
DELIVERY="$PROJECT/custom_research_report_20261023"
OUTDIR="$OUT/11_prepared_delivery_self_audit_v1"

[[ -s "$DELIVERY/P42_Vedita_18_WGS_1023.md" ]] || { echo "Prepared report is missing" >&2; exit 2; }
[[ ! -e "$OUTDIR" ]] || { echo "Refusing to overwrite existing stage: $OUTDIR" >&2; exit 2; }
"$PY" "$SCRIPT_DIR/11_validate_prepared_delivery.py" --project "$PROJECT" --delivery "$DELIVERY" --out-dir "$OUTDIR"
[[ -s "$OUTDIR/.complete" && -s "$OUTDIR/delivery_self_audit.json" ]] || { echo "Delivery self-audit incomplete" >&2; exit 2; }
