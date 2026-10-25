#!/usr/bin/env bash
# Revalidate the already built P42 client package after the Figure 12 correction.
set -euo pipefail

PROJECT=/home/gao/projects_2026Q4/42_Vedita_18_wgs_advanced_analysis
LOG_DIR="$PROJECT/logs/20261024_depth_interpretation_revision"
LOG="$LOG_DIR/final_validation_$(date +%Y%m%dT%H%M%S%z).log"
PYTHON=/Work_bio/gao/configs/.conda/envs/regular_bioinfo/bin/python

mkdir -p "$LOG_DIR"
exec >>"$LOG" 2>&1
printf 'START %s\n' "$(date --iso-8601=seconds)"
"$PYTHON" "$PROJECT/scripts/20261024_depth_interpretation_revision/03_validate_revised_client_package.py" \
  --source-package "$PROJECT/superseded_for_comparison/client_package_full/v3_20261024_remove_inapplicable_60x_gate/custom_research_report_20261024" \
  --delivery "$PROJECT/custom_research_report_20261024" \
  --audit-dir "$PROJECT/analysis/20261024_depth_interpretation_client_package_self_audit" \
  --samtools /Work_bio/gao/configs/.conda/envs/regular_bioinfo/bin/samtools \
  --bcftools /Work_bio/gao/configs/.conda/envs/regular_bioinfo/bin/bcftools
printf 'SUCCESS %s\n' "$(date --iso-8601=seconds)"
