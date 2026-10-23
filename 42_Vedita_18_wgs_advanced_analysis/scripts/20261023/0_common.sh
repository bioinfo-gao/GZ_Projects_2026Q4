#!/usr/bin/env bash
# Canonical constants for the 2026-10-23 full-cohort, repeat-aware P42 rerun.
# This file must be sourced, never executed as a pipeline itself.
set -euo pipefail

PROJECT=/home/gao/projects_2026Q4/42_Vedita_18_wgs_advanced_analysis
P41=/mnt/ex_8T_SSD/projects_2026M7/41_Vedict_18_wgs
SAMPLES="$PROJECT/scripts/0_sample_table.csv"
OUT="$PROJECT/analysis/20261023_fullcohort_repeat_resolve"
LOGDIR="$PROJECT/logs/20261023_fullcohort_repeat_resolve"
REF=/Work_bio/references/Homo_sapiens/GRCh38/GATK.GRCh38/Homo_sapiens_assembly38.fasta
PY=/Work_bio/gao/configs/.conda/envs/regular_bioinfo/bin/python
SAMTOOLS=/Work_bio/gao/miniforge3/bin/samtools
MOSDEPTH=/Work_bio/gao/configs/.conda/envs/regular_bioinfo/bin/mosdepth

TARGET_CHROM=chr11
TARGET_START=114604000
TARGET_END=114625500
DISCOVERY_START=114554000
DISCOVERY_END=114675500
TARGET_REGION="${TARGET_CHROM}:${TARGET_START}-${TARGET_END}"
DISCOVERY_REGION="${TARGET_CHROM}:${DISCOVERY_START}-${DISCOVERY_END}"

mkdir -p "$OUT" "$LOGDIR"

if [[ ! -s "$SAMPLES" ]]; then
  echo "Missing canonical sample table: $SAMPLES" >&2
  exit 2
fi
if [[ ! -s "$REF" || ! -s "${REF}.fai" ]]; then
  echo "Missing reference or fai: $REF" >&2
  exit 2
fi
if [[ $(($(wc -l < "$SAMPLES") - 1)) -ne 18 ]]; then
  echo "Expected exactly 18 data rows in $SAMPLES" >&2
  exit 2
fi

