#!/usr/bin/env bash
# P45 shared settings (adapted from P44 scripts/0_common.sh; snapshot in from_P44/). Source from every stage script.
set -euo pipefail

PROJ=/home/gao/projects_2026Q4/45_Vedita_18_wgs_combined_Advanced
P41=/home/gao/projects_2026M7/41_Vedict_18_wgs
P42=/home/gao/projects_2026Q4/42_Vedita_18_wgs_advanced_analysis
P44=/home/gao/projects_2026Q4/44_wgs_advanced_analysis_from_P41_by_cc
SAMPLE_TABLE=$PROJ/scripts/0_sample_table.csv
SCRIPTS=$PROJ/scripts/20261025

REF=/Work_bio/references/Homo_sapiens/GRCh38/GATK.GRCh38/Homo_sapiens_assembly38.fasta
BWA_IDX=/Work_bio/references/Homo_sapiens/GRCh38/GATK.GRCh38/BWAmem2Index/Homo_sapiens_assembly38.fasta

ENV=/Work_bio/gao/configs/.conda/envs/regular_bioinfo/bin
ASM_ENV=/Work_bio/gao/configs/.conda/envs/p44_assembly_20261024/bin
SV_ENV=/Work_bio/gao/configs/.conda/envs/p42_sv_validation_20261023/bin   # read-only use
SAMTOOLS=$ENV/samtools
PY=$ENV/python

# Target locus and extraction window: identical to P44 round 2 (W = target +/- 50 kb, user D3 2026-10-25)
TARGET_CHR=chr11; TARGET_START=114604000; TARGET_END=114625500
W_PAD=50000
W_START=$((TARGET_START - W_PAD)); W_END=$((TARGET_END + W_PAD))
W_REGION=${TARGET_CHR}:${W_START}-${W_END}

# Inherited (copied, read-only) P44 / P42 products, and P45's own outputs
P44_RUN=$P44/analysis/20261025_read_rescue          # stage1-8, W +/-50 kb
P44_SUPPORT=$P44/analysis/20261024_read_rescue/stage9_supporting   # W-independent
INH=$PROJ/analysis/20261025_inherited
OUT=$PROJ/analysis/20261025_p45
LOG=$PROJ/logs/20261025
mkdir -p "$INH" "$OUT" "$LOG"

samples_all() { awk -F, 'NR>1{print $2}' "$SAMPLE_TABLE"; }
samples_reported() { awk -F, 'NR>1 && $1<=10{print $2}' "$SAMPLE_TABLE"; }
bam_of() { awk -F, -v s="$1" 'NR>1 && $2==s{print $12}' "$SAMPLE_TABLE"; }
