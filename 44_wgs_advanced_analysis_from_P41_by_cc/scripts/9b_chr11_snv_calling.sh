#!/usr/bin/env bash
# Stage 9b-1 — joint germline SNV/indel calling on chr11 for S01-S10 (bcftools mpileup/call, multiallelic caller).
# Purpose: B-allele frequency / LOH along chr11 to interpret the S04/S05 chr11q depth step (supports Analysis 1).
# Scope per user decision Q3 (2026-10-24): chr11 only, not genome-wide.
# 5 Mb chunks, 8 in parallel (8 threads). Output: analysis/<RUN>/stage9_supporting/chr11_snv/
source "$(dirname "$0")/0_common.sh"
D=$OUT/stage9_supporting/chr11_snv; mkdir -p "$D/chunks"
BAMLIST=$D/bams_S01_S10.txt
awk -F, 'NR>1 && $1<=10{print $12}' "$SAMPLE_TABLE" > "$BAMLIST"
LEN=$(awk '$1=="chr11"{print $2}' "$REF.fai")
seq 0 5000000 "$LEN" | while read s; do e=$((s+5000000)); [[ $e -gt $LEN ]] && e=$LEN; echo "chr11:$((s+1))-$e"; done > "$D/chunks.txt"
call_chunk() {
  local r=$1 o="$D/chunks/${1//[:-]/_}.vcf.gz"
  [[ -s $o.csi ]] && return 0
  $ENV/bcftools mpileup -f "$REF" -r "$r" -b "$BAMLIST" -a FORMAT/AD,FORMAT/DP -q 20 -Q 20 --max-depth 1000 -Ou 2>/dev/null \
    | $ENV/bcftools call -mv -Oz -o "$o" && $ENV/bcftools index "$o"
}
export -f call_chunk; export D REF BAMLIST ENV
xargs -P 8 -I{} bash -c 'set -e; call_chunk {}' < "$D/chunks.txt"
ls "$D"/chunks/*.vcf.gz | sort -V > "$D/chunk_files.txt"
$ENV/bcftools concat -f "$D/chunk_files.txt" -Ou | $ENV/bcftools norm -f "$REF" -m -any -Oz -o "$D/chr11_S01_S10.raw.vcf.gz"
$ENV/bcftools index -t "$D/chr11_S01_S10.raw.vcf.gz"
$ENV/bcftools view -i 'QUAL>=30' -Oz -o "$D/chr11_S01_S10.QUAL30.vcf.gz" "$D/chr11_S01_S10.raw.vcf.gz"
$ENV/bcftools index -t "$D/chr11_S01_S10.QUAL30.vcf.gz"
echo "raw: $($ENV/bcftools index -n $D/chr11_S01_S10.raw.vcf.gz)  QUAL30: $($ENV/bcftools index -n $D/chr11_S01_S10.QUAL30.vcf.gz)"
touch "$D/.done"
