#!/usr/bin/env bash
# Stage G2a2 — parallel DELLY re-genotyping by site-chromosome groups (added 2026-10-25 15:55).
# Why: the single-process genotyping of S05 in G2a ran ~20x slower than S04 (1.2 vs 19 MB/s of BAM read; ETA ~15 h).
# Method: sites.bcf is split by the chromosome of each record (CHROM column; BND partner may be elsewhere), each part is
# genotyped with the identical delly sr -v command, parts are concatenated and sorted.
# Validation (must pass before the split result is used): the same split run on S04 must reproduce S04's single-run
# genotype BCF record-for-record (bcftools view -H, excluding nothing) — written to split_validation_S04.txt.
# Usage: G2a2_split_genotype.sh <sample> <outdir>   (does not touch the running G2a process or its outputs)
source "$(dirname "$0")/0_common.sh"
S=${1:?sample}; O=${2:?outdir}; mkdir -p "$O/parts"
D=$OUT/G2_dld1_genome_wide; EXCL=$P41/refs/delly_exclude_altdecoy_hla.tsv
CH=$($ENV/bcftools view -H "$D/sites.bcf" | cut -f1 | sort -u)
part() {
  local c=$1 o="$O/parts/$1"   # $1, not $c: both are expanded before local assigns
  [[ -s $o.geno.bcf ]] && return 0
  $ENV/bcftools view -r "$c" -O b -o "$o.sites.bcf" "$D/sites.bcf" && $ENV/bcftools index -f "$o.sites.bcf"
  OMP_NUM_THREADS=1 $ENV/delly sr -g "$REF" -x "$EXCL" -v "$o.sites.bcf" -o "$o.geno.tmp.bcf" "$(bam_of "$S")" > "$o.log" 2>&1 \
    || { echo "delly failed $S $c"; return 1; }
  mv "$o.geno.tmp.bcf" "$o.geno.bcf"; rm -f "$o.geno.tmp.bcf.csi"
}
export -f part bam_of; export O D REF EXCL ENV SAMPLE_TABLE S
echo "$CH" | xargs -P 14 -I{} bash -c 'set -e; part {}'
n=$(ls "$O"/parts/*.geno.bcf | wc -l); [[ $n -eq $(echo "$CH" | wc -l) ]] || { echo "parts incomplete $n"; exit 1; }
$ENV/bcftools concat -a -O u $(for c in $CH; do echo "$O/parts/$c.geno.bcf"; done) 2>/dev/null \
  | $ENV/bcftools sort -O b -o "$O/$S.geno.bcf" - && $ENV/bcftools index -f "$O/$S.geno.bcf"
echo "$(date '+%F %T') split genotype $S done: $($ENV/bcftools view -H "$O/$S.geno.bcf" | wc -l) records"
