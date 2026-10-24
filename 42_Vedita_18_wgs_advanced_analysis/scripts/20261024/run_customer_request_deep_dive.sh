#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
source "${SCRIPT_DIR}/6_customer_request_common.sh"

[[ ! -e "${OUT}" ]] || { echo "Refusing to overwrite existing analysis directory: ${OUT}" >&2; exit 1; }
mkdir -p "${LOGDIR}"
LOG_FILE="${LOGDIR}/customer_request_deep_dive_$(date +%Y%m%dT%H%M%S%z).log"
[[ ! -e "${LOG_FILE}" ]] || { echo "Refusing to overwrite existing log: ${LOG_FILE}" >&2; exit 1; }
exec >"${LOG_FILE}" 2>&1

echo "[$(date --iso-8601=seconds)] Starting P42 customer-request completion analysis"
"${PY}" "${SCRIPT_DIR}/6_customer_request_deep_dive.py" \
  --sample-table "${SAMPLES}" --p41 "${P41}" --summary "${SUMMARY}" --matrix "${MATRIX}" \
  --caller-input "${CALLER_INPUT}" --caller-records "${CALLER_RECORDS}" --rescue-ledger "${RESCUE}" \
  --raw-audit "${RAW_AUDIT}" --out-dir "${OUT}"
echo "[$(date --iso-8601=seconds)] Completed P42 customer-request completion analysis"
