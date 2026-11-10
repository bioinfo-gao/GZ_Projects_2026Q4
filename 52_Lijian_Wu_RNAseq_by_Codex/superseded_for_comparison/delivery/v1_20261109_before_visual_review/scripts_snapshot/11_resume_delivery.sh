#!/usr/bin/env bash
set -Eeuo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
P=/home/gao/projects_2026Q4/52_Lijian_Wu_RNAseq_by_Codex
trap 'echo "DELIVERY_RESUME_FAILED line=$LINENO $(date -Is)" >> "$P/logs/MASTER_STATUS"' ERR
python3 -u "$P/scripts/8a_render_figure_review.py"
python3 -u "$P/scripts/9_build_delivery.py" 2026-11-09
echo "DELIVERY_BUILT_AWAITING_VISUAL_VERIFICATION $(date -Is)" >> "$P/logs/MASTER_STATUS"
