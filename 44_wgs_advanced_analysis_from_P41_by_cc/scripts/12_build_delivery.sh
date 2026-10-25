#!/usr/bin/env bash
# Stage 12 — build the supplementary client delivery folder from canonical outputs (copy, never move).
# Usage: 12_build_delivery.sh <YYYY-MM-DD>   (date = generation date; re-run with the review-completion date
#        later: the existing folder is RENAMED, never duplicated — report skill rule).
# Hard checks (exit 1): exactly one delivery folder; exactly one .md in it; Report Date line equals the date;
# no @@REPORT_DATE@@ left; no CJK characters in any text deliverable; every file named in report section 9 exists.
source "$(dirname "$0")/0_common.sh"
DATE="${1:?usage: 12_build_delivery.sh YYYY-MM-DD}"; YMD="${DATE//-/}"; MMDD="${YMD:4:4}"
NAME=Chr11_Locus_Junction_Read_Rescue_Report
DELIV=$PROJ/custom_research_report_supplementary_$YMD
REPORT=${NAME}_${MMDD}.md
existing=$(find "$PROJ" -maxdepth 1 -type d -name 'custom_research_report_*' | head -1)
if [[ -n "$existing" && "$existing" != "$DELIV" ]]; then
  mv "$existing" "$DELIV"
  old=$(find "$DELIV" -maxdepth 1 -name "${NAME}_*.md" | head -1); [[ -n "$old" ]] && mv "$old" "$DELIV/$REPORT"
fi
mkdir -p "$DELIV"/{figures/junction_views,tables,supporting_analyses/{depth_boundary,chr11_BAF,chr11_SNV},genome_wide_SV_VCF,sequences,locus_BAM}
T=$OUT/stage10_tables; F=$OUT/stage11_figures; S9=$OUT/stage9_supporting
cp "$F"/P44_Fig[1-4]_*.{png,pdf} "$DELIV/figures/"
cp "$F"/junctions/P44_Fig5_*.{png,pdf} "$DELIV/figures/junction_views/"
cp "$T"/junction_table_target.tsv "$T"/junction_table_window_context.tsv "$T"/junction_by_sample_long.tsv \
   "$T"/rescue_source_summary.tsv "$DELIV/tables/"
cp "$OUT/stage4_clusters/filter_impact.tsv" "$DELIV/tables/evidence_filter_impact.tsv"
cp "$OUT/stage8_concordance/concordance_long.tsv" "$DELIV/tables/caller_concordance_by_sample.tsv"
cp "$S9"/depth_boundary/{boundary_estimates.tsv,depth_bins_200bp.tsv} "$DELIV/supporting_analyses/depth_boundary/"
cp "$S9"/chr11_baf/{baf_segments.tsv,retained_homolog_S04_vs_S05.tsv,baf_sites.tsv.gz} "$DELIV/supporting_analyses/chr11_BAF/"
cp "$S9"/chr11_snv/chr11_S01_S10.QUAL30.vcf.gz* "$DELIV/supporting_analyses/chr11_SNV/"
cp -r "$S9"/vcf_package/{manta,delly,md5sums.txt,vcf_inventory.tsv} "$DELIV/genome_wide_SV_VCF/"
cp "$OUT/stage6_genotype/templates.fa" "$DELIV/sequences/junction_and_reference_templates.fa"
cp "$OUT/stage5_assembly/all_contigs.fa" "$DELIV/sequences/local_assembly_contigs.fa"
for s in $(awk -F, 'NR>1 && $1<=10{print $2}' "$SAMPLE_TABLE"); do
  cp "$OUT/stage3_rescue/$s/rescue.bam" "$DELIV/locus_BAM/$s.locus_rescue.bam"
  cp "$OUT/stage3_rescue/$s/rescue.bam.bai" "$DELIV/locus_BAM/$s.locus_rescue.bam.bai"
done
sed "s/@@REPORT_DATE@@/$DATE/g" "$PROJ/docs/_client_report_source.md" > "$DELIV/$REPORT"
( cd "$DELIV" && find . -type f ! -name md5sums_delivery.txt -print0 | sort -z | xargs -0 md5sum > md5sums_delivery.txt )

# ---- hard checks
fail() { echo "CHECK FAILED: $*"; exit 1; }
[[ $(find "$PROJ" -maxdepth 1 -type d -name 'custom_research_report_*' | wc -l) -eq 1 ]] || fail "more than one delivery folder"
[[ $(find "$DELIV" -maxdepth 1 -name '*.md' | wc -l) -eq 1 ]] || fail "not exactly one report .md"
grep -q "^\*\*Report Date:\*\* $DATE\$" "$DELIV/$REPORT" || fail "Report Date line"
grep -q '@@REPORT_DATE@@' "$DELIV/$REPORT" && fail "placeholder left"
grep -rlP '[\x{4e00}-\x{9fff}]' "$DELIV" --include='*.md' --include='*.tsv' && fail "CJK text in deliverable"
for f in $(grep -oE '`[A-Za-z0-9_./-]+\.(tsv|png|pdf|fa|bam|vcf\.gz|gz|txt)`' "$DELIV/$REPORT" | tr -d '`' | sort -u); do
  [[ "$f" == *'<'* ]] && continue
  ls "$DELIV"/$f >/dev/null 2>&1 || ls "$DELIV"/*/$f >/dev/null 2>&1 || ls "$DELIV"/*/*/$f >/dev/null 2>&1 || fail "report names missing file: $f"
done
echo "delivery OK: $DELIV ($(find "$DELIV" -type f | wc -l) files, $(du -sh "$DELIV" | cut -f1))"
