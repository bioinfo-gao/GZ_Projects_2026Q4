#!/usr/bin/env bash
# P44 shared settings. Source from every stage script.
set -euo pipefail

PROJ=/home/gao/projects_2026Q4/44_wgs_advanced_analysis_from_P41_by_cc
P41=/home/gao/projects_2026M7/41_Vedict_18_wgs
SAMPLE_TABLE=$PROJ/scripts/0_sample_table.csv

REF=/Work_bio/references/Homo_sapiens/GRCh38/GATK.GRCh38/Homo_sapiens_assembly38.fasta
BWA_IDX=/Work_bio/references/Homo_sapiens/GRCh38/GATK.GRCh38/BWAmem2Index/Homo_sapiens_assembly38.fasta

ENV=/Work_bio/gao/configs/.conda/envs/regular_bioinfo/bin
ASM_ENV=/Work_bio/gao/configs/.conda/envs/p44_assembly_20261024/bin
SV_ENV=/Work_bio/gao/configs/.conda/envs/p42_sv_validation_20261023/bin   # read-only use (svaba)
SAMTOOLS=$ENV/samtools
PY=$ENV/python

# Target locus (client) and extraction window W = target +/- 5 kb (plan section 3)
TARGET_CHR=chr11; TARGET_START=114604000; TARGET_END=114625500
W_START=114599000; W_END=114630500
W_REGION=${TARGET_CHR}:${W_START}-${W_END}

# Run-dated output roots (never overwrite: each stage writes under its own dated dir)
RUN_DATE=${P44_RUN_DATE:-20261024}
OUT=$PROJ/analysis/${RUN_DATE}_read_rescue
LOG=$PROJ/logs/${RUN_DATE}_read_rescue
mkdir -p "$OUT" "$LOG"

# All 18 BAMs: S01-S10 reported, S11-S18 artefact background panel (user Q2, 2026-10-24)
samples_all() { awk -F, 'NR>1{print $2}' "$SAMPLE_TABLE"; }
bam_of() { awk -F, -v s="$1" 'NR>1 && $2==s{print $12}' "$SAMPLE_TABLE"; }
