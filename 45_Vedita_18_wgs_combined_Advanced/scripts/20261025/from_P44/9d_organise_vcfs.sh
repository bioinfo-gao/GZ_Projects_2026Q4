#!/usr/bin/env bash
# Stage 9d — organise the existing genome-wide P41 caller outputs for S01-S10 (user decision Q3, 2026-10-24:
# organise existing Manta/DELLY VCFs only; no new genome-wide calling). Read-only on P41; writes copies here.
#   manta/<S>.manta.diploidSV.unfiltered.vcf.gz   all diploidSV records (any FILTER)
#   manta/<S>.manta.diploidSV.PASS.vcf.gz         FILTER=PASS subset
#   manta/<S>.manta.candidateSV.vcf.gz            candidate-stage SVs (pre-scoring)
#   delly/<S>.delly.unfiltered.vcf.gz             delly sr discovery, all records
#   delly/<S>.delly.PASS.vcf.gz                   FILTER=PASS subset
# Each output is indexed; record counts written to vcf_inventory.tsv; md5 in md5sums.txt.
source "$(dirname "$0")/0_common.sh"
D=$OUT/stage9_supporting/vcf_package; mkdir -p "$D/manta" "$D/delly"
BC=$ENV/bcftools
printf "sample\tfile\trecords\n" > "$D/vcf_inventory.tsv"
for S in $(awk -F, 'NR>1 && $1<=10{print $2}' "$SAMPLE_TABLE"); do
  V=$P41/analysis/sv/manta_gw/$S/results/variants
  cp "$V/diploidSV.vcf.gz" "$D/manta/$S.manta.diploidSV.unfiltered.vcf.gz"
  cp "$V/candidateSV.vcf.gz" "$D/manta/$S.manta.candidateSV.vcf.gz"
  $BC view -f PASS -Oz -o "$D/manta/$S.manta.diploidSV.PASS.vcf.gz" "$V/diploidSV.vcf.gz"
  $BC view -Oz -o "$D/delly/$S.delly.unfiltered.vcf.gz" "$P41/analysis/sv/delly/discover/$S.delly.bcf"
  $BC view -f PASS -Oz -o "$D/delly/$S.delly.PASS.vcf.gz" "$P41/analysis/sv/delly/discover/$S.delly.bcf"
  for f in "$D"/manta/$S.*.vcf.gz "$D"/delly/$S.*.vcf.gz; do
    $BC index -f -t "$f"
    printf "%s\t%s\t%s\n" "$S" "$(basename "$f")" "$($BC index -n "$f")" >> "$D/vcf_inventory.tsv"
  done
done
( cd "$D" && md5sum manta/*.vcf.gz delly/*.vcf.gz > md5sums.txt )
touch "$D/.done"
echo "stage9d OK"
