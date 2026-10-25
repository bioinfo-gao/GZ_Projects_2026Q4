#!/usr/bin/env bash
# Stage 0 — wait for P44 round 2 (W +/-50 kb) CHAIN_DONE, then copy (not link) the products P45 reuses unchanged.
# Read-only on P44/P42. Re-entrant: files whose md5 already matches are skipped. Output: analysis/20261025_inherited/
source "$(dirname "$0")/0_common.sh"
P44_CHAINLOG=$P44/logs/20261025_read_rescue/chain.log
until grep -q CHAIN_DONE "$P44_CHAINLOG" 2>/dev/null; do
  if grep -qiE "incomplete|Traceback" "$P44_CHAINLOG" 2>/dev/null; then echo "P44 chain reports failure:"; tail -5 "$P44_CHAINLOG"; exit 2; fi
  sleep 120
  # chain runs under set -e: a failed stage ends the tmux session without CHAIN_DONE
  if ! tmux has-session -t p44r2_chain 2>/dev/null && ! grep -q CHAIN_DONE "$P44_CHAINLOG"; then
    echo "P44 tmux p44r2_chain ended without CHAIN_DONE:"; tail -5 "$P44_CHAINLOG"; exit 2; fi
done
echo "$(date '+%F %T') P44 round-2 CHAIN_DONE seen"

# 1. provenance: P44 scripts at import time vs P45 snapshot (from_P44/) taken 2026-10-25 12:05
( cd "$P44/scripts" && md5sum *.py *.sh ) | sort -k2 > "$INH/p44_scripts_md5_at_import.txt"
( cd "$SCRIPTS/from_P44" && md5sum *.py *.sh ) | sort -k2 > "$INH/p44_scripts_md5_snapshot.txt"
diff "$INH/p44_scripts_md5_snapshot.txt" "$INH/p44_scripts_md5_at_import.txt" > "$INH/p44_scripts_diff.txt" \
  && echo "P44 scripts unchanged since snapshot" || echo "WARNING: P44 scripts changed since snapshot (see p44_scripts_diff.txt)"

# 2. copy list: <source path> <destination relative to INH>
L=$INH/copy_list.tsv; : > "$L"
for st in stage1_inventory stage2_scan stage3_rescue stage4_clusters stage5_assembly stage6_genotype stage7_annotation stage8_concordance; do
  [[ -d $P44_RUN/$st ]] || { echo "missing $P44_RUN/$st"; exit 1; }
  echo -e "$P44_RUN/$st\tP44_r2/$st" >> "$L"
done
echo -e "$P44_SUPPORT\tP44_r1_stage9_supporting" >> "$L"
echo -e "$P42/analysis/20261023_fullcohort_repeat_resolve/07_external_sv_dld1_normalization_v1\tP42/07_external_sv_dld1_normalization_v1" >> "$L"
for f in fixed_candidate_summary.tsv fixed_candidate_dld1_matrix.tsv c2_read_end_junction_resolution.tsv c2_junction_resolution_summary.tsv \
         candidate_manta_delly_concordance.tsv read_rescue_candidate_impact.tsv clone_specific_comparison_with_criteria.tsv \
         primer_specificity.tsv assay_design_selection.tsv standard_evidence_autosomal_depth.tsv customer_request_coverage.tsv \
         data_inventory.tsv external_input_status.tsv raw_support_alignment_audit.tsv; do
  echo -e "$P42/custom_research_report_20261024/results/$f\tP42/results/$f" >> "$L"
done
for f in virtual_junction_templates.fasta bam_header_sanitization_audit.tsv vcf_header_sanitization_audit.tsv delly_conversion_audit.tsv; do
  echo -e "$P42/custom_research_report_20261024/methods/$f\tP42/methods/$f" >> "$L"
done

# 3. copy (cp -a: preserves mtimes, never hard links) + manifest
M=$INH/import_manifest.tsv
echo -e "source\tdest\tbytes\tmd5_source\tmd5_dest" > "$M.tmp"
while IFS=$'\t' read -r src dst; do
  mkdir -p "$INH/$(dirname "$dst")"
  if [[ -d $src ]]; then cp -a --no-preserve=links "$src/." "$INH/$dst/"; else cp -a "$src" "$INH/$dst"; fi
  while IFS= read -r -d '' f; do
    rel=${f#$src}; d="$INH/$dst$rel"
    ms=$(md5sum < "$f" | cut -c1-32); md=$(md5sum < "$d" | cut -c1-32)
    echo -e "$f\t${d#$INH/}\t$(stat -c %s "$f")\t$ms\t$md" >> "$M.tmp"
  done < <(if [[ -d $src ]]; then find "$src" -type f -print0; else printf '%s\0' "$src"; fi)
done < "$L"
mv "$M.tmp" "$M"
bad=$(awk -F'\t' 'NR>1 && $4!=$5' "$M" | wc -l)
echo "files: $(($(wc -l < "$M")-1))  md5 mismatches: $bad  total: $(du -sh "$INH" | cut -f1)"
[[ $bad -eq 0 ]] || exit 1
touch "$INH/.done"; echo "$(date '+%F %T') STAGE0 DONE"
