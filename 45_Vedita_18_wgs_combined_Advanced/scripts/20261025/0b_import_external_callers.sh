#!/usr/bin/env bash
# Stage 0b — copy P42's GRIDSS (18-sample joint) and SvABA (DLD1 S01-S05 multisample, v5) VCFs for Stage 8x.
# Both were run by P42 on chr11:114,554,000-114,675,500 (= P45 W; P42 regions/discovery_window_0based.bed).
# SvABA v5 is the version P42 used in 07_external_sv_dld1_normalization_v1 (audit.json input_sample_sets).
# Read-only on P42; independent of P44, so it runs immediately. Output: analysis/20261025_inherited/P42/external_callers/
source "$(dirname "$0")/0_common.sh"
X=$P42/analysis/20261023_fullcohort_repeat_resolve
D=$INH/P42/external_callers; mkdir -p "$D"
SRC=(
  "$X/06_external_sv_validation_v3/gridss_joint/gridss_joint_18samples.vcf.gz"
  "$X/06_external_sv_validation_v3/gridss_joint/gridss_joint_18samples.vcf.gz.tbi"
  "$X/06_external_sv_validation_v3/regions/discovery_window_0based.bed"
  "$X/06_external_svaba_dld1_v5_multisample/svaba_dld1/svaba_dld1_targeted_multisample.svaba.germline.sv.vcf.gz"
  "$X/06_external_svaba_dld1_v5_multisample/svaba_dld1/svaba_dld1_targeted_multisample.svaba.somatic.sv.vcf.gz"
  "$X/06_external_svaba_dld1_v5_multisample/svaba_dld1/svaba_dld1_targeted_multisample.svaba.unfiltered.germline.sv.vcf.gz"
  "$X/06_external_svaba_dld1_v5_multisample/svaba_dld1/svaba_dld1_targeted_multisample.svaba.unfiltered.somatic.sv.vcf.gz"
  "$X/06_external_svaba_dld1_v5_multisample/svaba_dld1/input_manifest.tsv"
)
echo -e "source\tdest\tbytes\tmd5_source\tmd5_dest" > "$D/import_manifest.tsv"
for f in "${SRC[@]}"; do
  cp -a "$f" "$D/"; d=$D/$(basename "$f")
  echo -e "$f\t${d#$INH/}\t$(stat -c %s "$f")\t$(md5sum < "$f" | cut -c1-32)\t$(md5sum < "$d" | cut -c1-32)" >> "$D/import_manifest.tsv"
done
awk -F'\t' 'NR>1 && $4!=$5{bad++} END{print "copied", NR-1, "md5 mismatches", bad+0; exit bad>0}' "$D/import_manifest.tsv"
touch "$D/.done"
