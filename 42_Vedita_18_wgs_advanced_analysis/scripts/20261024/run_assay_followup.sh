#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
source "${SCRIPT_DIR}/0_common.sh"

if [[ -e "${OUT}" ]]; then
    echo "Refusing to reuse existing analysis directory: ${OUT}" >&2
    exit 1
fi
if [[ -e "${DELIVERY}" ]]; then
    echo "Refusing to overwrite existing delivery directory: ${DELIVERY}" >&2
    exit 1
fi

mkdir -p "${LOGDIR}"
LOG_FILE="${LOGDIR}/p42_assay_followup_$(date +%Y%m%dT%H%M%S%z).log"
if [[ -e "${LOG_FILE}" ]]; then
    echo "Refusing to overwrite existing run log: ${LOG_FILE}" >&2
    exit 1
fi
exec >"${LOG_FILE}" 2>&1

echo "[$(date --iso-8601=seconds)] Starting P42 P41-Section-7.3 assay follow-up"
echo "Output: ${OUT}"
echo "Delivery: ${DELIVERY}"

"${PRIMER_PY}" "${SCRIPT_DIR}/1_design_targeted_assays.py" \
    --summary "${SUMMARY}" \
    --ledger "${LEDGER}" \
    --normalized "${NORMALIZED}" \
    --samples "${SAMPLES}" \
    --reference "${REF}" \
    --samtools "${SAMTOOLS}" \
    --out-dir "${OUT}"

"${FIGURE_PY}" "${SCRIPT_DIR}/2_check_primer_specificity.py" \
    --unchecked-pairs "${OUT}/primer_pairs_unchecked.tsv" \
    --template-audit "${OUT}/candidate_template_audit.tsv" \
    --reference "${REF}" \
    --makeblastdb "${MAKEBLASTDB}" \
    --blastn "${BLASTN}" \
    --out-dir "${OUT}"

mkdir -p "${OUT}/figures"
"${FIGURE_PY}" "${SCRIPT_DIR}/3_make_assay_followup_figure.py" \
    --selection "${OUT}/assay_design_selection.tsv" \
    --out-prefix "${OUT}/figures/P42_Fig13_P41_7_3_assay_followup"

"${FIGURE_PY}" "${SCRIPT_DIR}/4_build_supplementary_delivery.py" \
    --project "${PROJECT}" \
    --out-dir "${OUT}" \
    --delivery "${DELIVERY}"

"${FIGURE_PY}" "${SCRIPT_DIR}/5_validate_supplementary_delivery.py" \
    --out-dir "${OUT}" \
    --delivery "${DELIVERY}" \
    --audit-dir "${OUT}/05_delivery_self_audit"

touch "${OUT}/.pipeline_complete"
echo "[$(date --iso-8601=seconds)] Completed P42 P41-Section-7.3 assay follow-up"
