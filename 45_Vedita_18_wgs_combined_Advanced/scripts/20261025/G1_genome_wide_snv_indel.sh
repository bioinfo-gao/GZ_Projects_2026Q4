#!/usr/bin/env bash
# Stage G1 — genome-wide joint germline SNV/indel calling, S01-S10 (user Q3 = b, 2026-10-25).
# Same caller and parameters as the P44 chr11 call (from_P44/9b_chr11_snv_calling.sh) so chr11 results stay comparable:
#   bcftools mpileup -q 20 -Q 20 --max-depth 1000 -a AD,DP | bcftools call -mv ; norm -m -any ; QUAL>=30 filtered set.
# Scope: primary assembly chr1-22, X, Y, M (alt/decoy/HLA/EBV contigs not called; stated in Methods).
# 5 Mb chunks, 20 in parallel (P44 chain may use 28 threads concurrently; total stays <= 56).
# Resumable: finished chunks (index present) are skipped. Output: analysis/20261025_p45/G1_snv_indel/
source "$(dirname "$0")/0_common.sh"
D=$OUT/G1_snv_indel; mkdir -p "$D/chunks"
BAMLIST=$D/bams_S01_S10.txt
awk -F, 'NR>1 && $1<=10{print $12}' "$SAMPLE_TABLE" > "$BAMLIST"
[[ $(wc -l < "$BAMLIST") -eq 10 ]] || { echo "expected 10 BAMs"; exit 1; }
CHROMS="$(printf 'chr%s ' {1..22} X Y M)"
for c in $CHROMS; do
  L=$(awk -v c="$c" '$1==c{print $2}' "$REF.fai"); [[ -n $L ]] || { echo "missing $c in fai"; exit 1; }
  for ((s=0; s<L; s+=5000000)); do e=$((s+5000000)); ((e>L)) && e=$L; echo "$c:$((s+1))-$e"; done
done > "$D/chunks.txt"
echo "$(date '+%F %T') chunks: $(wc -l < "$D/chunks.txt")"
call_chunk() {
  local r=$1 o="$D/chunks/${1//[:-]/_}.vcf.gz"
  [[ -s $o.csi ]] && return 0
  $ENV/bcftools mpileup -f "$REF" -r "$r" -b "$BAMLIST" -a FORMAT/AD,FORMAT/DP -q 20 -Q 20 --max-depth 1000 -Ou 2> "$o.mpileup.log" \
    | $ENV/bcftools call -mv -Oz -o "$o" && $ENV/bcftools index "$o"
}
export -f call_chunk; export D REF BAMLIST ENV
xargs -P 20 -I{} bash -c 'set -e; call_chunk {}' < "$D/chunks.txt"
n_exp=$(wc -l < "$D/chunks.txt"); n_ok=$(ls "$D"/chunks/*.vcf.gz.csi | wc -l)
[[ $n_ok -eq $n_exp ]] || { echo "chunks incomplete: $n_ok/$n_exp"; exit 1; }
# concat in chunks.txt order (genomic order), not lexical
while read r; do echo "$D/chunks/${r//[:-]/_}.vcf.gz"; done < "$D/chunks.txt" > "$D/chunk_files.txt"
$ENV/bcftools concat -f "$D/chunk_files.txt" -Ou | $ENV/bcftools norm -f "$REF" -m -any -Oz -o "$D/S01_S10.genome.raw.vcf.gz"
$ENV/bcftools index -t "$D/S01_S10.genome.raw.vcf.gz"
$ENV/bcftools view -i 'QUAL>=30' -Oz -o "$D/S01_S10.genome.QUAL30.vcf.gz" "$D/S01_S10.genome.raw.vcf.gz"
$ENV/bcftools index -t "$D/S01_S10.genome.QUAL30.vcf.gz"
$ENV/bcftools stats "$D/S01_S10.genome.QUAL30.vcf.gz" > "$D/S01_S10.genome.QUAL30.stats.txt"
# consistency check vs the P44 chr11 call (same parameters => chr11 record counts should match exactly)
c11_new=$($ENV/bcftools view -H -r chr11 "$D/S01_S10.genome.QUAL30.vcf.gz" | wc -l)
c11_p44=$($ENV/bcftools view -H "$P44_SUPPORT/chr11_snv/chr11_S01_S10.QUAL30.vcf.gz" | wc -l)
echo -e "metric\tvalue\nraw_records\t$($ENV/bcftools index -n $D/S01_S10.genome.raw.vcf.gz)\nQUAL30_records\t$($ENV/bcftools index -n $D/S01_S10.genome.QUAL30.vcf.gz)\nchr11_QUAL30_G1\t$c11_new\nchr11_QUAL30_P44\t$c11_p44" | tee "$D/summary.tsv"
touch "$D/.done"
echo "$(date '+%F %T') G1 DONE"
