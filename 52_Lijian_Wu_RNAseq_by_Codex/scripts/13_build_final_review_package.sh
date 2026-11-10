#!/usr/bin/env bash
set -Eeuo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
P=/home/gao/projects_2026Q4/52_Lijian_Wu_RNAseq_by_Codex
trap 'echo "FINAL_PACKAGE_FAILED line=$LINENO $(date -Is)" >> "$P/logs/MASTER_STATUS"' ERR
python3 -u "$P/scripts/6b_add_volcano_zoom.py"
python3 -u "$P/scripts/9_build_delivery.py" 2026-11-09
/usr/bin/google-chrome --headless --disable-gpu --disable-background-networking --disable-sync --no-first-run --user-data-dir=/tmp/p52_chrome_review_20261109 --window-size=1600,1200 --timeout=20000 --screenshot="$P/logs/multiqc_browser_20261109.png" "file://$P/custom_research_report_20261109/QC/multiqc/star_salmon/multiqc_report.html" > "$P/logs/multiqc_browser_20261109.log" 2>&1
echo "FINAL_PACKAGE_AUDITED $(date -Is)" >> "$P/logs/MASTER_STATUS"
