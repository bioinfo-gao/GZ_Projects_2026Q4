#!/usr/bin/env bash
# Constants for the client-question completion supplement; all inputs are read-only.
set -euo pipefail

PROJECT=/home/gao/projects_2026Q4/42_Vedita_18_wgs_advanced_analysis
P41=/mnt/ex_8T_SSD/projects_2026M7/41_Vedict_18_wgs
PY=/Work_bio/gao/configs/.conda/envs/regular_bioinfo/bin/python
SAMPLES="$PROJECT/scripts/0_sample_table.csv"
PRIOR="$PROJECT/analysis/20261023_fullcohort_repeat_resolve"
ASSAY="$PROJECT/analysis/20261024_p41_7_3_assay_followup_v4_coordinate_specificity"
OUT="$PROJECT/analysis/20261024_customer_request_response_v1"
LOGDIR="$PROJECT/logs/20261024_customer_request_response"

SUMMARY="$PRIOR/09_fixed_candidate_regenotyping_v1/fixed_candidate_summary.tsv"
MATRIX="$PRIOR/09_fixed_candidate_regenotyping_v1/fixed_candidate_dld1_matrix.tsv"
CALLER_INPUT="$PRIOR/02_caller_audit/caller_input_audit.tsv"
CALLER_RECORDS="$PRIOR/02_caller_audit/caller_records_in_discovery_window.tsv"
RESCUE="$PRIOR/05_local_assembly/inputs_v1/unmapped_mate_rescue/unmapped_mate_recovery_ledger.tsv"
RAW_AUDIT="$ASSAY/raw_support_alignment_audit.tsv"

for f in "$SAMPLES" "$SUMMARY" "$MATRIX" "$CALLER_INPUT" "$CALLER_RECORDS" "$RESCUE" "$RAW_AUDIT"; do
    [[ -s "$f" ]] || { echo "Missing required input: $f" >&2; exit 2; }
done
