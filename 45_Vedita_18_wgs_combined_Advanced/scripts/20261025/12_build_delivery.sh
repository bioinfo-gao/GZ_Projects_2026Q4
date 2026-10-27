#!/usr/bin/env bash
# Stage 12 (P45; from P44 12_build_delivery.sh) — build the single complete client package from canonical outputs (copy only).
# Usage: 12_build_delivery.sh <YYYY-MM-DD>
# P45 changes: one complete package custom_research_report_<YMD>/ (not supplementary); internal junction IDs (P44J...)
# rewritten to client IDs (J...) in every copied table/FASTA; GRIDSS/SvABA locus VCFs, genome-wide SNV/indel VCFs (G1) and
# the genome-wide DLD1 comparison (G2) added; P42-only rows (comparison primer pair) not delivered; the pre-sanitisation
# md5sums.txt of the VCF package is not copied (the delivery md5 covers every file); chr11-only SNV VCF replaced by the
# genome-wide set. A previous build of the same date is MOVED to superseded_for_comparison/delivery_prebuilds/, never deleted.
# Hard checks (exit 1): one report .md; Report Date line; no placeholder; no CJK; no internal project numbers in any
# text deliverable; every file named in the report exists.
source "$(dirname "$0")/0_common.sh"
DATE="${1:?usage: 12_build_delivery.sh YYYY-MM-DD}"; YMD="${DATE//-/}"; MMDD="${YMD:4:4}"
NAME=Chr11_Locus_and_DLD1_Clone_Analysis_Report
DELIV=$PROJ/custom_research_report_$YMD
REPORT=${NAME}_${MMDD}.md
if [[ -e $DELIV ]]; then
  mkdir -p "$PROJ/superseded_for_comparison/delivery_prebuilds"
  mv "$DELIV" "$PROJ/superseded_for_comparison/delivery_prebuilds/$(basename "$DELIV")_moved_$(date +%Y%m%d_%H%M%S)"
fi
mkdir -p "$DELIV"/{methods,figures/junction_views,tables,sequences,locus_BAM,supporting_analyses/{depth_boundary,chr11_BAF,SV_verification,DLD1_sub_lineage}} \
         "$DELIV"/{genome_wide_SV_VCF/locus_GRIDSS_SvABA,genome_wide_SNV_indel_VCF,genome_wide_DLD1_comparison}
T=$OUT/stage10_tables; F=$OUT/stage11_figures; R2=$INH/P44_r2; S9=$INH/P44_r1_stage9_supporting
G1=$OUT/G1_snv_indel; G2=$OUT/G2_dld1_genome_wide; X=$INH/P42/external_callers
cid() { sed -E 's/\bP44J([0-9]+)/J\1/g'; }
# 1027: figures are numbered in the order they appear in the report (Results), not in script order. Analysis outputs keep
# their script names; the report number is assigned here. FIGMAP is the single place where the mapping lives.
declare -A FIGMAP=([Fig1_locus_overview]=Fig1_locus_overview [Fig3_rescue_sources]=Fig2_rescue_sources
                   [Fig0_key_events_summary]=Fig3_key_events_summary [Fig4_before_after_and_classes]=Fig5_before_after_and_classes
                   [Fig2_chr11_BAF]=Fig6_chr11_BAF)
for src in "${!FIGMAP[@]}"; do for ext in png pdf; do cp "$F/$src.$ext" "$DELIV/figures/${FIGMAP[$src]}.$ext"; done; done
for p in "$F"/junctions/Fig5_*.{png,pdf}; do b=$(basename "$p"); cp "$p" "$DELIV/figures/junction_views/Fig4_${b#Fig5_}"; done
$PY - "$F/figure_manifest.json" "$DELIV/figures/figure_manifest.json" \
    $(for k in "${!FIGMAP[@]}"; do printf '%s=%s ' "$k" "${FIGMAP[$k]}"; done) <<'PYEOF'
import json, sys
m = dict(a.split("=", 1) for a in sys.argv[3:])          # FIGMAP from the shell, so the mapping exists only once
man = json.load(open(sys.argv[1]))
for x in man:
    f = x["file"].replace("junctions/Fig5_", "junction_views/Fig4_")
    stem, ext = f.rsplit(".", 1)
    x["file"] = f"{m.get(stem, stem)}.{ext}"
json.dump(man, open(sys.argv[2], "w"), indent=1)
PYEOF
cp "$T"/junction_table_target.tsv "$T"/junction_table_window_context.tsv "$T"/junction_by_sample_long.tsv \
   "$T"/rescue_source_summary.tsv "$T"/key_junction_read_evidence.tsv "$DELIV/tables/"
cid < "$R2/stage4_clusters/filter_impact.tsv" > "$DELIV/tables/evidence_filter_impact.tsv"
cid < "$OUT/stage8x_concordance/concordance_long_4callers.tsv" > "$DELIV/tables/caller_concordance_by_sample.tsv"
cid < "$OUT/stage9x_coherence/coherence_tests.tsv" > "$DELIV/tables/structural_coherence_tests.tsv"
cp "$OUT"/stage14_primers/{primer_design.tsv,primer_candidates.tsv} "$DELIV/tables/"
awk -F'\t' -v OFS='\t' 'NR==1{$1="pair"} NR>1{$1=($1=="P45_selected")?"selected":"alternative_both_primers_in_flanks"} 1' \
  "$OUT/stage14_primers/primer_binding_sites.tsv" > "$DELIV/tables/primer_binding_sites.tsv"
cid < "$R2/stage6_genotype/templates.fa" > "$DELIV/sequences/junction_and_reference_templates.fa"
cid < "$R2/stage5_assembly/all_contigs.fa" > "$DELIV/sequences/local_assembly_contigs.fa"
for s in $(samples_reported); do
  cp "$R2/stage3_rescue/$s/rescue.bam" "$DELIV/locus_BAM/$s.locus_rescue.bam"
  cp "$R2/stage3_rescue/$s/rescue.bam.bai" "$DELIV/locus_BAM/$s.locus_rescue.bam.bai"
done
cp "$S9"/depth_boundary/{boundary_estimates.tsv,depth_bins_200bp.tsv} "$DELIV/supporting_analyses/depth_boundary/"
cp "$S9"/chr11_baf/{baf_segments.tsv,retained_homolog_S04_vs_S05.tsv,baf_sites.tsv.gz} "$DELIV/supporting_analyses/chr11_BAF/"
# 2026-10-27 checks cited in Sections 7.7, 7.9 and 8.3 (scripts/20261027/15_*, 16_*)
cp "$PROJ"/analysis/20261027_sv_verification/*.tsv "$DELIV/supporting_analyses/SV_verification/"
cp "$PROJ"/analysis/20261027_dld1_lineage/*.tsv "$DELIV/supporting_analyses/DLD1_sub_lineage/"
cp -r "$S9"/vcf_package/{manta,delly,vcf_inventory.tsv} "$DELIV/genome_wide_SV_VCF/"
cp "$X"/gridss_joint_18samples.vcf.gz "$DELIV/genome_wide_SV_VCF/locus_GRIDSS_SvABA/GRIDSS_joint_18samples.window.vcf.gz"
for k in germline somatic; do
  cp "$X/svaba_dld1_targeted_multisample.svaba.$k.sv.vcf.gz" "$DELIV/genome_wide_SV_VCF/locus_GRIDSS_SvABA/SvABA_DLD1_S01_S05.window.$k.filtered.vcf.gz"
  cp "$X/svaba_dld1_targeted_multisample.svaba.unfiltered.$k.sv.vcf.gz" "$DELIV/genome_wide_SV_VCF/locus_GRIDSS_SvABA/SvABA_DLD1_S01_S05.window.$k.unfiltered.vcf.gz"
done
cp "$G1"/S01_S10.genome.{raw,QUAL30}.vcf.gz "$G1/S01_S10.genome.QUAL30.stats.txt" "$DELIV/genome_wide_SNV_indel_VCF/"
cp "$G2"/{dld1_genome_wide_events.tsv,dld1_dox_specific_events.tsv,summary_counts.tsv,somatic_samples.tsv} "$DELIV/genome_wide_DLD1_comparison/"
$ENV/bcftools view -Oz -o "$DELIV/genome_wide_DLD1_comparison/dld1_S01_S05.delly_genotyped.vcf.gz" "$G2/dld1_S01_S05.somatic_tagged.bcf"
sed "s/@@REPORT_DATE@@/$DATE/g" "$PROJ/docs/_client_report_source.md" > "$DELIV/$REPORT"
$PY "$SCRIPTS/12a_sanitize_headers.py" "$DELIV" || { echo "header sanitisation failed"; exit 1; }   # also indexes VCFs
$PY "$SCRIPTS/12c_delly_conversion_audit.py" "$DELIV" || { echo "DELLY conversion audit failed"; exit 1; }
for v in $(find "$DELIV" -name '*.vcf.gz'); do [[ -e $v.tbi || -e $v.csi ]] || $ENV/bcftools index -f -t "$v"; done
$PY "$SCRIPTS/12b_delivery_tables.py" "$DELIV" || { echo "delivery tables failed"; exit 1; }
( cd "$DELIV" && find . -type f ! -name md5sums_delivery.txt -print0 | sort -z | xargs -0 md5sum > md5sums_delivery.txt )

# ---- hard checks
fail() { echo "CHECK FAILED: $*"; exit 1; }
[[ $(find "$DELIV" -maxdepth 1 -name '*.md' | wc -l) -eq 1 ]] || fail "not exactly one report .md"
grep -q "^\*\*Report Date:\*\* $DATE\$" "$DELIV/$REPORT" || fail "Report Date line"
grep -q '@@REPORT_DATE@@' "$DELIV/$REPORT" && fail "placeholder left"
grep -rlP '[\x{4e00}-\x{9fff}]' "$DELIV" --include='*.md' --include='*.tsv' --include='*.json' && fail "CJK text in deliverable"
grep -rlE '\bP4[0-9]J?[0-9]*\b|P4[0-9]_' "$DELIV" --include='*.md' --include='*.tsv' --include='*.json' --include='*.fa' && fail "internal project number in deliverable"
for f in $(grep -oE '`[A-Za-z0-9_./-]+\.(tsv|png|pdf|fa|bam|vcf\.gz|gz|txt|json)`' "$DELIV/$REPORT" | tr -d '`' | sort -u); do
  [[ "$f" == *'<'* ]] && continue
  ls "$DELIV"/$f >/dev/null 2>&1 || ls "$DELIV"/*/$f >/dev/null 2>&1 || ls "$DELIV"/*/*/$f >/dev/null 2>&1 || fail "report names missing file: $f"
done
echo "delivery OK: $DELIV ($(find "$DELIV" -type f | wc -l) files, $(du -sh "$DELIV" | cut -f1))"
