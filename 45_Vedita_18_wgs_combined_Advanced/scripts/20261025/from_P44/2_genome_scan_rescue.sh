#!/usr/bin/env bash
# Stage 2 — one streaming pass over each whole BAM (read-only) to find records that belong to the locus
# but are not stored inside W:
#   (b) mate (RNEXT/PNEXT) inside W, (c) any SA segment inside W, (d) any XA alternative hit inside W.
# Also writes (a) all records stored inside W and the both-mates-unmapped pairs for source (e).
# The regex is a superset of W (exact interval re-checked in Stage 3); nothing is dropped here except
# records that match none of the conditions.
# Usage: 2_genome_scan_rescue.sh <sample_id>      (runner: 2_run_scan_all.sh)
source "$(dirname "$0")/0_common.sh"
S=$1; B=$(bam_of "$S"); D=$OUT/stage2_scan/$S; mkdir -p "$D"
[[ -f $D/.done ]] && { echo "[$S] already done"; exit 0; }

if [[ $W_PAD -eq 5000 ]]; then   # round 1 expression (kept verbatim for reproducibility)
  RX='chr11,[+-]?(114599[0-9]{3}|1146[0-2][0-9]{4}|114630[0-4][0-9]{2}|1146305[0-9]{2})'
else                              # superset 114,500,000-114,699,999; exact W re-checked in Stage 3
  [[ $W_START -ge 114500000 && $W_END -le 114699999 ]] || { echo "regex superset does not cover W"; exit 1; }
  RX='chr11,[+-]?114[56][0-9]{5}'
fi
EXPR="(rnext==\"chr11\" && pnext>=${W_START} && pnext<=${W_END}) || (rnext==\"=\" && rname==\"chr11\" && pnext>=${W_START} && pnext<=${W_END}) || [SA]=~\"$RX\" || [XA]=~\"$RX\""

echo "[$S] $(date '+%F %T') start; BAM=$B"
$SAMTOOLS view -@ 3 -b -o "$D/window.bam" "$B" "$W_REGION"
$SAMTOOLS view -@ 3 -b -e "$EXPR" -o "$D/scan.bam" "$B"
$SAMTOOLS view -@ 3 -b -f 12 -o "$D/unmapped_pairs.bam" "$B" '*'
for f in window scan unmapped_pairs; do
  printf "%s\t%s\t%s\n" "$S" "$f" "$($SAMTOOLS view -c "$D/$f.bam")"
done > "$D/counts.tsv"
cat "$D/counts.tsv"
echo "$EXPR" > "$D/filter_expression.txt"
touch "$D/.done"
echo "[$S] $(date '+%F %T') done"
