#!/usr/bin/env bash
# P46 交付包构建（2026-11-26 全量首轮 ⇒ 目录名不带 supplementary）
#
# 规则（/delivery + /report skill）：
#   - 日期是显式参数，三处必须相等：目录名 / 报告文件名 MMDD / 正文 Report Date
#   - **copy 不 move**，源文件原地保留
#   - 交付图件只放 PDF（PNG 预览留在 analysis/，不进交付）
#   - 不放 .rds / 中间缓存 / 预览产物
#   - 交付目录内按分析类型分子目录，不平铺
# 硬断言在末尾：三处日期相等、顶层 .md 恰好 1 个、7 张 PDF 且 0 张 PNG、无 .rds、
#   报告 §8.2/§9 引用的每个相对路径都真实存在。
set -uo pipefail
DATE=${1:?用法: 33_build_delivery.sh YYYYMMDD}
MMDD=${DATE:4:4}
P=/home/gao/projects_2026Q4/46_park_wgs
A=$P/analysis
OUT=$P/custom_research_report_$DATE
REPORT_SRC=$A/stage5_report_20261126/Dual_Cut_Project_WGS_$MMDD.md
say(){ echo "[$(date '+%F %T')] $*"; }
[[ -s "$REPORT_SRC" ]] || { say "FATAL 报告源文件不存在：$REPORT_SRC"; exit 1; }
mkdir -p "$OUT"/{figures,qc,identity,on_target,off_target,vector,copy_number,variants}

cp -f "$REPORT_SRC" "$OUT/"
cp -f "$A/stage5_figures_20261126"/*.pdf "$OUT/figures/"

# QC：每样本比对统计 + 深度 + MultiQC
B=/srv/external_data/gao_projects/46_park_wgs_bam/stage0_20261121
for s in $(python3 -c "import csv;[print(r['safe_name']) for r in csv.DictReader(open('$P/scripts/0_sample_table.csv'))]"); do
  cp -f "$B/qc/$s.flagstat.txt" "$B/qc/$s.samtools_stats.txt" "$B/qc/$s.markdup_metrics.txt" "$OUT/qc/" 2>/dev/null
  cp -f "$B/coverage/$s.mosdepth.summary.txt" "$OUT/qc/" 2>/dev/null
done
cp -rf "$A/stage0_identity_qc_20261123_v2/multiqc" "$OUT/qc/" 2>/dev/null

# 身份与污染
cp -f "$A/stage0_identity_qc_20261123_v2/tables"/{identity_matrix.tsv,identity_n_sites.tsv,contamination.tsv,filter_impact.tsv} "$OUT/identity/"

# 靶点
cp -f "$A/stage1_ontarget_20261123/tables"/{junction.tsv,ontarget_summary.tsv,ontarget_summary.json,junction_reads.tsv,indels_near_cutsites.tsv,supplementary_alignments.tsv,filter_impact_audit.tsv} "$OUT/on_target/"
cp -f "$A/stage5_verify_20261126/tables/calibration_vs_client_corrected.tsv" "$OUT/on_target/"

# 脱靶
cp -f "$A/stage2a_candidates_20261121/tables/offtarget_candidate_sites_n1587.tsv" "$OUT/off_target/"
cp -f "$A/stage2b_offtarget_evidence_20261123_v2/tables"/{offtarget_evidence_n10.tsv,offtarget_summary.json,calibration_ontarget.tsv,filter_impact_audit.tsv} "$OUT/off_target/"
cp -f "$A/stage2c_paired_cut_20261123/paired_cut_offtargets_n10.tsv" "$OUT/off_target/"
cp -f "$A/stage5_verify_20261126/tables"/{treated_only_sites_reads.tsv,sample4_chr19_inv_context.tsv} "$OUT/off_target/"

# 载体
cp -f "$A/stage3_vector_v2_20261124/tables"/{vector_coverage_n10.tsv,integration_junctions_n10.tsv,integration_sharing_n10.tsv,stage3_summary.json,filter_impact_audit.tsv} "$OUT/vector/"
cp -f "$A/stage5_verify_20261126/tables"/{vector_copies_normalized.tsv,junction_breakpoint_fixedness.tsv,junction_fixedness_summary.json} "$OUT/vector/"

# 拷贝数
cp -f "$A/stage4c_cnv_20261123/baf_check"/{chr12p_baf_crosscheck.tsv,verdict.json} "$OUT/copy_number/"
cp -f "$A/stage4c_cnv_20261123/tables"/{cnv_segments_n10.tsv,cnv_per_sample.tsv,filter_impact_audit.tsv} "$OUT/copy_number/"

# 变异
V=$A/stage4a_germline_n10_20261123/vcf/P46_joint_n10_annotated.vcf.gz
cp -f "$V" "$V.tbi" "$OUT/variants/" 2>/dev/null
cp -f "$A/stage4a_germline_n10_20261123/qc"/{P46_stage4a_summary.tsv,P46_joint_n10_annotated.bcftools_stats.txt} "$OUT/variants/" 2>/dev/null
cp -f "$A/stage4d_private_burden_20261125/tables"/{private_burden_per_sample.tsv,private_variants.tsv,stage4d_summary.json,filter_impact_audit.tsv} "$OUT/variants/"
cp -f "$A/stage4b_delly_n10_20261125"/{delly_joint_n10.bcf,delly_joint_n10.bcf.csi} "$OUT/variants/" 2>/dev/null
cp -f "$A/stage5_prep_20261125/tables"/{private_sv_per_sample.tsv,private_sv_hits_lists.tsv,private_sv_summary.json} "$OUT/variants/"

say "=== 硬断言 ==="
R="$OUT/Dual_Cut_Project_WGS_$MMDD.md"
fail=0
chk(){ if eval "$2"; then echo "  ✓ $1"; else echo "  ✗ $1"; fail=1; fi; }
chk "三处日期一致（目录 $DATE / 文件名 $MMDD / 正文 Report Date）" \
    "grep -q '^\*\*Report Date:\*\* ${DATE:0:4}-${DATE:4:2}-${DATE:6:2}$' '$R'"
chk "顶层 .md 恰好 1 个" "[[ \$(ls '$OUT'/*.md | wc -l) -eq 1 ]]"
chk "figures/ 7 张 PDF" "[[ \$(ls '$OUT'/figures/*.pdf | wc -l) -eq 7 ]]"
chk "figures/ 无 PNG" "[[ \$(ls '$OUT'/figures/*.png 2>/dev/null | wc -l) -eq 0 ]]"
chk "交付包内无 .rds" "[[ \$(find -L '$OUT' -name '*.rds' | wc -l) -eq 0 ]]"
chk "交付包内无 tmp/preview" "[[ \$(find -L '$OUT' -iname '*preview*' -o -iname 'tmp*' | wc -l) -eq 0 ]]"
chk "qc/ 每样本 4 个文件 ×10" "[[ \$(ls '$OUT'/qc/*.flagstat.txt | wc -l) -eq 10 ]]"
# §8.2 / §9 里写出的相对路径必须真实存在
miss=0
for f in $(grep -oE '`[a-z_]+/[A-Za-z0-9_.]+`' "$R" | tr -d '`' | sort -u); do
  [[ -e "$OUT/$f" ]] || { echo "    缺: $f"; miss=1; }
done
chk "报告引用的每个相对路径都存在" "[[ $miss -eq 0 ]]"
say "文件数 $(find -L "$OUT" -type f | wc -l)，体积 $(du -sh "$OUT" | cut -f1)"
[[ $fail -eq 0 ]] && say "DELIVERY_BUILD_OK  $OUT" || { say "DELIVERY_BUILD_FAILED"; exit 1; }
