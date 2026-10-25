#!/usr/bin/env bash
# Rebuild the single P42 client package after correcting the depth interpretation.
# The old generated-but-unreviewed package remains immutable in superseded_for_comparison/.
set -euo pipefail

PROJECT=/home/gao/projects_2026Q4/42_Vedita_18_wgs_advanced_analysis
REVISION_DIR="$PROJECT/analysis/20261024_depth_interpretation_revision_v1"
SOURCE="$PROJECT/superseded_for_comparison/client_package_full/v3_20261024_remove_inapplicable_60x_gate/custom_research_report_20261024"
DESTINATION="$PROJECT/custom_research_report_20261024"
AUDIT_DIR="$PROJECT/analysis/20261024_depth_interpretation_client_package_self_audit"
LOG_DIR="$PROJECT/logs/20261024_depth_interpretation_revision"
PYTHON=/Work_bio/gao/configs/.conda/envs/regular_bioinfo/bin/python
SAMTOOLS=/Work_bio/gao/configs/.conda/envs/regular_bioinfo/bin/samtools
BCFTOOLS=/Work_bio/gao/configs/.conda/envs/regular_bioinfo/bin/bcftools

mkdir -p "$LOG_DIR"
LOG="$LOG_DIR/rebuild_and_validate_$(date +%Y%m%dT%H%M%S%z).log"
exec >>"$LOG" 2>&1

printf 'START %s\n' "$(date --iso-8601=seconds)"
printf 'SOURCE=%s\nDESTINATION=%s\nAUDIT_DIR=%s\n' "$SOURCE" "$DESTINATION" "$AUDIT_DIR"

"$PYTHON" "$PROJECT/scripts/20261024_depth_interpretation_revision/02_build_revised_client_package.py" \
  --source-package "$SOURCE" \
  --destination "$DESTINATION" \
  --depth-dir "$REVISION_DIR" \
  --decision-table "$REVISION_DIR/analysis_decision_table.tsv"

"$PYTHON" "$PROJECT/scripts/20261024_depth_interpretation_revision/03_validate_revised_client_package.py" \
  --source-package "$SOURCE" \
  --delivery "$DESTINATION" \
  --audit-dir "$AUDIT_DIR" \
  --samtools "$SAMTOOLS" \
  --bcftools "$BCFTOOLS"

printf 'SUCCESS %s\n' "$(date --iso-8601=seconds)"
