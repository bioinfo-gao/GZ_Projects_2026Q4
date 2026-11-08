#!/usr/bin/env bash
# =============================================================================
# P51 — build / re-stamp the client delivery folder (adapted from proj30 9_build_delivery.sh).
#
# DATES (house rule 2026-08-09): first build = today; while review is running DO NOT re-date;
# once the user says the review is complete, re-run with that date: folder name, report filename
# and "Report Date:" are re-stamped together, by RENAMING the one existing folder (never a copy).
# No `rm` anywhere: a report carrying an older date is MOVED to superseded_for_comparison/report/.
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

# ---- older-dated reports in the source folder -> archive (one live report only) ----
for f in "${SRC}"/Crtc2_Line_Treatment_RNAseq_Report_*.md; do
  [ -e "$f" ] || continue
  [ "$(basename "$f")" = "$REPORT_NAME" ] && continue
  mkdir -p "$ARCH"; mv "$f" "${ARCH}/$(basename "$f" .md)_from_Data_Analysis_moved$(date +%Y%m%d%H%M).md"
  echo "ARCHIVED older source report: $(basename "$f")"
done

# ---- the date has ONE source: regenerate the report with it ----
DELIVERY_DATE="$ISO" "$RSCRIPT" "${SCRIPT_DIR}/6_generate_research_report.R" > "${PROJ_DIR}/logs/report_${YMD}.log" 2>&1 \
  || { echo "ABORT: report generation failed (logs/report_${YMD}.log)"; tail -20 "${PROJ_DIR}/logs/report_${YMD}.log"; exit 1; }
[ -s "${SRC}/${REPORT_NAME}" ] || { echo "ABORT: report not produced"; exit 1; }

# ---- audit the source before copying ----
"$RSCRIPT" "${SCRIPT_DIR}/8_pre_delivery_audit.R" > "${PROJ_DIR}/logs/audit_source_${YMD}.log" 2>&1 \
  || { echo "ABORT: source audit failed"; grep -v PASS "${PROJ_DIR}/logs/audit_source_${YMD}.log" | tail -20; exit 1; }

# ---- one delivery folder only: re-stamp in place ----
mapfile -t EXISTING < <(ls -d "${PROJ_DIR}"/custom_research_report_[0-9]* 2>/dev/null)
[ "${#EXISTING[@]}" -le 1 ] || { echo "ABORT: more than one delivery folder:"; printf '  %s\n' "${EXISTING[@]}"; exit 1; }
if [ "${#EXISTING[@]}" -eq 1 ] && [ "${EXISTING[0]}" != "$DEST" ]; then
  echo "RE-STAMP: $(basename "${EXISTING[0]}") -> $(basename "$DEST") (rename in place)"
  mv "${EXISTING[0]}" "$DEST" || { echo "ABORT: rename failed"; exit 1; }
elif [ "${#EXISTING[@]}" -eq 0 ]; then
  echo "FIRST BUILD: $(basename "$DEST")"
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
