#!/usr/bin/env bash
# Stage G2a — genome-wide DLD1 clone comparison (client Analysis 2; user Q4 = b, 2026-10-25), DELLY germline->somatic workflow:
#   1. union of the P41 per-sample DELLY discovery sites of S01-S05 (delly merge; quality/AF gates relaxed to -y 0 -a 0.05
#      so the site universe is not pre-filtered: every site is re-genotyped and judged by its per-sample read counts)
#   2. re-genotype every site in each of the 5 BAMs (delly sr -v; same exclude list as P41: alt/decoy/HLA/chrM)
#   3. bcftools merge -> one 5-sample BCF
#   4. delly filter -f somatic with S04,S05 = tumor, S01-S03 = control, -t (TAG, never remove: filter-impact stays visible)
# Per-event labels with explicit criteria are made in G2b. Threads: 3 samples x OMP 2 (G1 and P44 run concurrently).
source "$(dirname "$0")/0_common.sh"
D=$OUT/G2_dld1_genome_wide; mkdir -p "$D/geno"
EXCL=$P41/refs/delly_exclude_altdecoy_hla.tsv
DISC=$P41/analysis/sv/delly/discover
DLD1=(S01_DLD1_pop S02_DLD1_c1_ctrl S03_DLD1_c2_ctrl S04_DLD1_c3_dox S05_DLD1_c4_dox)
die() { echo "FAILED: $*"; exit 1; }
if [[ ! -s $D/sites.bcf ]]; then
  $ENV/delly merge -y 0 -a 0.05 -v 5 -o "$D/sites.tmp.bcf" $(for s in "${DLD1[@]}"; do echo "$DISC/$s.delly.bcf"; done) \
    > "$D/merge.log" 2>&1 || die "delly merge"
  [[ -s $D/sites.tmp.bcf ]] || die "sites.bcf empty"
  mv "$D/sites.tmp.bcf" "$D/sites.bcf"; mv "$D/sites.tmp.bcf.csi" "$D/sites.bcf.csi" 2>/dev/null || $ENV/bcftools index -f "$D/sites.bcf"
fi
echo "$(date '+%F %T') sites: $($ENV/bcftools view -H "$D/sites.bcf" | wc -l)"
geno_one() {
  local s=$1 o="$D/geno/$1.geno.bcf"
  [[ -s $o ]] && return 0
  OMP_NUM_THREADS=2 $ENV/delly sr -g "$REF" -x "$EXCL" -v "$D/sites.bcf" -o "$o.tmp.bcf" "$(bam_of "$s")" > "$D/geno/$s.log" 2>&1 \
    || { echo "delly genotype failed: $s"; return 1; }
  [[ -s $o.tmp.bcf ]] || { echo "empty genotype BCF: $s"; return 1; }
  mv "$o.tmp.bcf" "$o"; mv "$o.tmp.bcf.csi" "$o.csi" 2>/dev/null || $ENV/bcftools index -f "$o"
  echo "$(date '+%F %T') genotyped $s"
}
export -f geno_one bam_of; export D REF EXCL ENV SAMPLE_TABLE
printf '%s\n' "${DLD1[@]}" | xargs -P 3 -I{} bash -c 'set -e; geno_one {}'
[[ $(ls "$D"/geno/*.geno.bcf | wc -l) -eq 5 ]] || die "genotype BCFs != 5"
$ENV/bcftools merge -m id -O b -o "$D/dld1_S01_S05.geno.bcf" $(for s in "${DLD1[@]}"; do echo "$D/geno/$s.geno.bcf"; done)
$ENV/bcftools index -f "$D/dld1_S01_S05.geno.bcf"
printf '%s\ttumor\n' S04_DLD1_c3_dox S05_DLD1_c4_dox > "$D/somatic_samples.tsv"
printf '%s\tcontrol\n' S01_DLD1_pop S02_DLD1_c1_ctrl S03_DLD1_c2_ctrl >> "$D/somatic_samples.tsv"
$ENV/delly filter -f somatic -t -s "$D/somatic_samples.tsv" -o "$D/dld1_S01_S05.somatic_tagged.bcf" "$D/dld1_S01_S05.geno.bcf" \
  > "$D/filter.log" 2>&1 || die "delly filter"
$ENV/bcftools index -f "$D/dld1_S01_S05.somatic_tagged.bcf"
echo -e "metric\tvalue\nsites_merged\t$($ENV/bcftools view -H $D/sites.bcf | wc -l)\ngenotyped_records\t$($ENV/bcftools view -H $D/dld1_S01_S05.geno.bcf | wc -l)\nsomatic_PASS\t$($ENV/bcftools view -H -f PASS $D/dld1_S01_S05.somatic_tagged.bcf | wc -l)" | tee "$D/summary.tsv"
touch "$D/.G2a_done"; echo "$(date '+%F %T') G2a DONE"
