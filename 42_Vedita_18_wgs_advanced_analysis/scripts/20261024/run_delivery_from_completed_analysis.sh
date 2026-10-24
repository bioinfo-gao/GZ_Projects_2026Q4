#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
source "${SCRIPT_DIR}/0_common.sh"

for required in \
    "${OUT}/.design_complete" \
    "${OUT}/.specificity_complete" \
    "${OUT}/assay_design_selection.tsv" \
    "${OUT}/figures/P42_Fig13_P41_7_3_assay_followup.png" \
    "${OUT}/figures/P42_Fig13_P41_7_3_assay_followup.pdf"; do
    [[ -f "${required}" ]] || { echo "Missing completed-analysis prerequisite: ${required}" >&2; exit 1; }
done
[[ ! -e "${DELIVERY}" ]] || { echo "Refusing to overwrite existing delivery directory: ${DELIVERY}" >&2; exit 1; }

mkdir -p "${LOGDIR}"
LOG_FILE="${LOGDIR}/p42_delivery_rebuild_$(date +%Y%m%dT%H%M%S%z).log"
[[ ! -e "${LOG_FILE}" ]] || { echo "Refusing to overwrite existing delivery-stage log: ${LOG_FILE}" >&2; exit 1; }
exec >"${LOG_FILE}" 2>&1

echo "[$(date --iso-8601=seconds)] Starting non-overwriting delivery build from completed v2 analysis"
"${FIGURE_PY}" "${SCRIPT_DIR}/4_build_supplementary_delivery.py" \
    --project "${PROJECT}" \
    --out-dir "${OUT}" \
    --delivery "${DELIVERY}"
"${FIGURE_PY}" "${SCRIPT_DIR}/5_validate_supplementary_delivery.py" \
    --out-dir "${OUT}" \
    --delivery "${DELIVERY}" \
    --audit-dir "${OUT}/05_delivery_self_audit_v2"
touch "${OUT}/.pipeline_complete"
echo "[$(date --iso-8601=seconds)] Completed supplementary delivery build and self-audit"
