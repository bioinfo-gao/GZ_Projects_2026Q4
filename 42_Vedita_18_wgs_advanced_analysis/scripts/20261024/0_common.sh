#!/usr/bin/env bash
# Constants for the 2026-10-24 incremental P41 §7.3 assay-design follow-up.
set -euo pipefail

PROJECT=/home/gao/projects_2026Q4/42_Vedita_18_wgs_advanced_analysis
P41=/mnt/ex_8T_SSD/projects_2026M7/41_Vedict_18_wgs
REF=/Work_bio/references/Homo_sapiens/GRCh38/GATK.GRCh38/Homo_sapiens_assembly38.fasta
SAMTOOLS=/Work_bio/gao/miniforge3/bin/samtools
MAKEBLASTDB=/Work_bio/gao/configs/.conda/envs/regular_bioinfo/bin/makeblastdb
BLASTN=/Work_bio/gao/configs/.conda/envs/regular_bioinfo/bin/blastn
PRIMER_PY=/Work_bio/gao/configs/.conda/envs/p42_assay_design_20261024/bin/python
FIGURE_PY=/Work_bio/gao/configs/.conda/envs/regular_bioinfo/bin/python
SAMPLES="$PROJECT/scripts/0_sample_table.csv"
S05_BAM=/mnt/ex_8T_SSD/projects_2026M7/41_Vedict_18_wgs/analysis/bam/S05_DLD1_c4_dox.markdup.bam
PRIOR_OUT="$PROJECT/analysis/20261023_fullcohort_repeat_resolve"
OUT="$PROJECT/analysis/20261024_p41_7_3_assay_followup_v4_coordinate_specificity"
LOGDIR="$PROJECT/logs/20261024_p41_7_3_assay_followup"
DELIVERY="$PROJECT/custom_research_report_supplementary_20261024_v5"

SUMMARY="$PRIOR_OUT/09_fixed_candidate_regenotyping_v1/fixed_candidate_summary.tsv"
LEDGER="$PRIOR_OUT/09_fixed_candidate_regenotyping_v1/raw_qname_evidence_ledger.tsv"
NORMALIZED="$PRIOR_OUT/07_external_sv_dld1_normalization_v1/normalized_external_sv_events.tsv"

for f in "$REF" "$REF.fai" "$SAMPLES" "$S05_BAM" "$S05_BAM.bai" "$SUMMARY" "$LEDGER" "$NORMALIZED"; do
  [[ -s "$f" ]] || { echo "Missing required input: $f" >&2; exit 2; }
done
