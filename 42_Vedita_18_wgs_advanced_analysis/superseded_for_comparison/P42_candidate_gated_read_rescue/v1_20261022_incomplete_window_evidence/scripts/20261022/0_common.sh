#!/usr/bin/env bash
set -euo pipefail

PROJ=/home/gao/projects_2026Q4/42_Vedita_18_wgs_advanced_analysis
RUN_DATE=20261022
RUN="$PROJ/analysis/$RUN_DATE"
SCRIPTS="$PROJ/scripts/$RUN_DATE"
SAMPLES="$PROJ/scripts/0_sample_table.csv"
P41=/mnt/ex_8T_SSD/projects_2026M7/41_Vedict_18_wgs
P41_EVIDENCE="$P41/analysis/region/evidence/target"
REF=/Work_bio/references/Homo_sapiens/GRCh38/GATK.GRCh38/Homo_sapiens_assembly38.fasta
BWA_REF=/Work_bio/references/Homo_sapiens/GRCh38/GATK.GRCh38/BWAmem2Index/Homo_sapiens_assembly38.fasta
ENVBIN=/Work_bio/gao/configs/.conda/envs/regular_bioinfo/bin
PY="$ENVBIN/python"
SAMTOOLS="$ENVBIN/samtools"
BWA="$ENVBIN/bwa-mem2"
TARGET=chr11:114604000-114625500
WINDOW=chr11:114554000-114675500
mkdir -p "$RUN" "$PROJ/logs"

die(){ echo "[$(date '+%F %T')] FATAL: $*" >&2; exit 1; }
need(){ [[ -s "$1" ]] || die "missing or empty: $1"; }
