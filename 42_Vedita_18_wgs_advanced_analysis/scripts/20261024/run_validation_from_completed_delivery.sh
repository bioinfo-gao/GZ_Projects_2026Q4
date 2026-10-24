#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
source "${SCRIPT_DIR}/0_common.sh"

[[ -f "${OUT}/.design_complete" ]] || { echo "Missing completed-analysis marker" >&2; exit 1; }
[[ -d "${DELIVERY}" ]] || { echo "Missing completed delivery directory" >&2; exit 1; }
AUDIT_DIR="${OUT}/05_delivery_self_audit_v3"
[[ ! -e "${AUDIT_DIR}" ]] || { echo "Refusing to overwrite existing validation audit: ${AUDIT_DIR}" >&2; exit 1; }

mkdir -p "${LOGDIR}"
LOG_FILE="${LOGDIR}/p42_delivery_validation_$(date +%Y%m%dT%H%M%S%z).log"
[[ ! -e "${LOG_FILE}" ]] || { echo "Refusing to overwrite existing validation log: ${LOG_FILE}" >&2; exit 1; }
exec >"${LOG_FILE}" 2>&1

echo "[$(date --iso-8601=seconds)] Starting non-overwriting supplementary delivery validation"
"${FIGURE_PY}" "${SCRIPT_DIR}/5_validate_supplementary_delivery.py" \
    --out-dir "${OUT}" \
    --delivery "${DELIVERY}" \
    --audit-dir "${AUDIT_DIR}"
touch "${OUT}/.pipeline_complete"
echo "[$(date --iso-8601=seconds)] Completed supplementary delivery validation"
