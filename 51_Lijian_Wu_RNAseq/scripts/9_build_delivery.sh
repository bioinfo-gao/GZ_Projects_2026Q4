#!/usr/bin/env bash
# =============================================================================
# P51 — build / re-stamp the client delivery folder (adapted from proj30 9_build_delivery.sh).
#
# VERSIONS (user directive 2026-10-27 in P45, reaffirmed 2026-11-10): each new version gets its OWN
# folder custom_research_report_<YYYYMMDD>/ and its own dated report; earlier folders and earlier-dated
# reports in Data_Analysis_* are never moved, renamed or overwritten, so versions compare at a glance.
# The folder of the version the user finally signs off is the submission folder. No `rm` anywhere.
#
# Usage: bash scripts/9_build_delivery.sh [YYYY-MM-DD]
# =============================================================================
set -uo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; PROJ_DIR="$(dirname "$SCRIPT_DIR")"
RSCRIPT=/Work_bio/gao/configs/.conda/envs/DE_R45/bin/Rscript
REVIEW_DATE="${1:-$(date +%F)}"
date -d "$REVIEW_DATE" >/dev/null 2>&1 || { echo "ABORT: bad date '$REVIEW_DATE'"; exit 1; }
YMD=$(date -d "$REVIEW_DATE" +%Y%m%d); MMDD=$(date -d "$REVIEW_DATE" +%m%d); ISO=$(date -d "$REVIEW_DATE" +%F)
DEST="${PROJ_DIR}/custom_research_report_${YMD}"
REPORT_NAME="Crtc2_Line_Treatment_RNAseq_Report_${MMDD}.md"
ARCH="${PROJ_DIR}/superseded_for_comparison/report"
mkdir -p "${PROJ_DIR}/logs"

mapfile -t SRCS < <(ls -d "${PROJ_DIR}"/Data_Analysis_[0-9]* 2>/dev/null)
[ "${#SRCS[@]}" -eq 1 ] || { echo "ABORT: need exactly one Data_Analysis_* (found ${#SRCS[@]})"; exit 1; }
SRC="${SRCS[0]}"
echo "Source: $(basename "$SRC") | Delivery date: ${ISO} | Target: $(basename "$DEST")"

# ---- earlier-dated reports in Data_Analysis_* stay where they are (one per version) ----
[ ! -e "${SRC}/${REPORT_NAME}" ] || [ -d "$DEST" ] || { echo "ABORT: ${REPORT_NAME} already exists but ${DEST##*/} does not; use a new date"; exit 1; }

# ---- the date has ONE source: regenerate the report with it ----
DELIVERY_DATE="$ISO" "$RSCRIPT" "${SCRIPT_DIR}/6_generate_research_report.R" > "${PROJ_DIR}/logs/report_${YMD}.log" 2>&1 \
  || { echo "ABORT: report generation failed (logs/report_${YMD}.log)"; tail -20 "${PROJ_DIR}/logs/report_${YMD}.log"; exit 1; }
[ -s "${SRC}/${REPORT_NAME}" ] || { echo "ABORT: report not produced"; exit 1; }

# ---- new version = new folder; earlier custom_research_report_* folders are never touched ----
# (Data_Analysis_* may hold several dated reports, so the audit runs on the delivery folder only, below.)
if [ ! -d "$DEST" ]; then
  echo "NEW VERSION FOLDER: $(basename "$DEST")"
  mkdir -p "$DEST"
  for d in QC Reads DE_PCA_Results Enrichment; do cp -r "${SRC}/${d}" "${DEST}/${d}"; done
  cp "${SRC}"/mouse_Gene_annotation_*.xlsx "${DEST}/"
else
  echo "REFRESH: folder already dated ${YMD}"
fi
for f in "${DEST}"/*.md; do
  [ -e "$f" ] || continue; [ "$(basename "$f")" = "$REPORT_NAME" ] && continue
  mkdir -p "$ARCH"; mv "$f" "${ARCH}/$(basename "$f" .md)_from_delivery_moved$(date +%Y%m%d%H%M).md"
  echo "ARCHIVED older delivered report: $(basename "$f")"
done
cp "${SRC}/${REPORT_NAME}" "${DEST}/${REPORT_NAME}"
for f in Crtc2_line_treatment_response_retention.csv CREB_CRTC_target_and_inflammatory_marker_panel.csv; do
  cp "${SRC}/DE_PCA_Results/${f}" "${DEST}/DE_PCA_Results/${f}"     # report-step outputs: keep in sync
done

# ---- checks on the delivery folder itself ----
FAIL=0
BODY_DATE=$(grep -m1 '^\*\*Report Date:\*\*' "${DEST}/${REPORT_NAME}" | grep -oE '[0-9]{4}-[0-9]{2}-[0-9]{2}')
[ "$(basename "$DEST")" = "custom_research_report_${YMD}" ] || { echo "FAIL folder name"; FAIL=1; }
[ "$BODY_DATE" = "$ISO" ] || { echo "FAIL Report Date ${BODY_DATE}"; FAIL=1; }
[ "$(ls "${DEST}"/*.md | wc -l)" -eq 1 ] || { echo "FAIL: not exactly one report"; FAIL=1; }
cmp -s "${SRC}/${REPORT_NAME}" "${DEST}/${REPORT_NAME}" || { echo "FAIL: delivered report != source"; FAIL=1; }
diff -rq "${SRC}/DE_PCA_Results" "${DEST}/DE_PCA_Results" >/dev/null || { echo "FAIL: DE_PCA_Results differ from source"; FAIL=1; }
diff -rq "${SRC}/Enrichment" "${DEST}/Enrichment" >/dev/null || { echo "FAIL: Enrichment differs from source"; FAIL=1; }
[ "$FAIL" -eq 0 ] || exit 1
"$RSCRIPT" "${SCRIPT_DIR}/8_pre_delivery_audit.R" "$DEST" > "${PROJ_DIR}/logs/audit_delivery_${YMD}.log" 2>&1 \
  || { echo "FAIL: delivery audit"; grep -v PASS "${PROJ_DIR}/logs/audit_delivery_${YMD}.log" | tail -20; exit 1; }
python3 /home/gao/.claude/hooks/checks/check_report.py "${DEST}/${REPORT_NAME}" || { echo "FAIL: check_report.py"; exit 1; }
echo "OK: dates = ${ISO}; one report; delivery == source; audit $(tail -2 "${PROJ_DIR}/logs/audit_delivery_${YMD}.log" | head -1)"
echo "Delivery: ${DEST}  Size: $(du -sh "$DEST" | cut -f1)  Files: $(find "$DEST" -type f | wc -l)"
