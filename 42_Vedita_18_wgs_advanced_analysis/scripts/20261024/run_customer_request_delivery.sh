#!/usr/bin/env bash
# Render, package, and self-audit only the new P42 customer-question outputs.
set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
source "${SCRIPT_DIR}/6_customer_request_common.sh"

DELIVERY="${PROJECT}/custom_research_report_supplementary_20261024_v8"
AUDIT_DIR="${OUT}/14_customer_request_delivery_self_audit_v8"
FIG_PREFIX="${OUT}/figures/P42_Fig14_Client_Request_Completion"
LOG_FILE="${LOGDIR}/customer_request_delivery_$(date +%Y%m%dT%H%M%S%z).log"

[[ ! -e "${DELIVERY}" ]] || { echo "Refusing to overwrite existing delivery: ${DELIVERY}" >&2; exit 1; }
[[ ! -e "${AUDIT_DIR}" ]] || { echo "Refusing to overwrite existing audit directory: ${AUDIT_DIR}" >&2; exit 1; }
[[ ! -e "${LOG_FILE}" ]] || { echo "Refusing to overwrite existing log: ${LOG_FILE}" >&2; exit 1; }
mkdir -p "${LOGDIR}"
exec >"${LOG_FILE}" 2>&1

echo "[$(date --iso-8601=seconds)] Starting P42 customer-question delivery generation"
if [[ -e "${FIG_PREFIX}.png" || -e "${FIG_PREFIX}.pdf" ]]; then
    [[ -s "${FIG_PREFIX}.png" && -s "${FIG_PREFIX}.pdf" ]] || { echo "Figure 14 exists only partially; refusing to overwrite it" >&2; exit 1; }
    echo "[$(date --iso-8601=seconds)] Reusing the already-rendered, geometry-checked Figure 14 source without overwriting it"
else
    "${PY}" "${SCRIPT_DIR}/7_make_customer_request_figure.py" --input-dir "${OUT}" --out-prefix "${FIG_PREFIX}"
fi
"${PY}" "${SCRIPT_DIR}/8_build_customer_request_delivery.py" --out-dir "${OUT}" --delivery "${DELIVERY}"
"${PY}" "${SCRIPT_DIR}/9_validate_customer_request_delivery.py" --out-dir "${OUT}" --delivery "${DELIVERY}" --audit-dir "${AUDIT_DIR}"
echo "[$(date --iso-8601=seconds)] Completed P42 customer-question delivery generation"
