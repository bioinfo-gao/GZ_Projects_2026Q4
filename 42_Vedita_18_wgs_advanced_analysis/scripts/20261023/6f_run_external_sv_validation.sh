#!/usr/bin/env bash
# P42 external validation: targeted GRIDSS2 and SvABA with preserved raw context.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/0_common.sh"

VALID="${P42_VALID_STAGE:-$OUT/06_external_sv_validation}"
# The validated GRIDSS reference is a completed shared prerequisite, not a
# per-attempt output; versioned validation attempts consume it read-only.
GRIDSS_REF="$OUT/06_external_sv_validation/gridss_reference_v1/Homo_sapiens_assembly38.fasta"
GRIDSS_ENV=/Work_bio/gao/configs/.conda/envs/p42_sv_validation_20261023
GRIDSS="$GRIDSS_ENV/bin/gridss"
EXTRACT="$GRIDSS_ENV/bin/gridss_extract_overlapping_fragments"
SVABA="$GRIDSS_ENV/bin/svaba"
AUDITOR="$SCRIPT_DIR/6f_audit_gridss_targeted_bam.py"
BED="$VALID/regions/discovery_window_0based.bed"
RAW="$VALID/gridss_targeted_inputs/raw_extractor"
CALLER="$VALID/gridss_targeted_inputs/caller_compatible"
AUDIT="$VALID/gridss_targeted_inputs/audit"
RUNLOG="$LOGDIR/06_external_sv_validation.run.log"

mkdir -p "$VALID/regions" "$RAW" "$CALLER" "$AUDIT" "$VALID/gridss_joint" "$VALID/svaba_dld1"
export PATH="$GRIDSS_ENV/bin:$PATH"
[[ -s "$GRIDSS_REF" && -s "$GRIDSS_REF.bwt" && -s "$GRIDSS_REF.sa" ]] || {
  echo "GRIDSS private reference/index incomplete: $GRIDSS_REF" >&2; exit 2;
}
[[ -x "$GRIDSS" && -x "$EXTRACT" && -x "$SVABA" && -s "$AUDITOR" ]] || {
  echo "Validation tool or audit script missing" >&2; exit 2;
}

# BED is 0-based half-open; it denotes the already registered 1-based discovery window.
printf 'chr11\t114553999\t114675500\tP42_discovery_window\n' > "$BED"

declare -a GRIDSS_BAMS=()
declare -a DLD1_CASE_BAMS=()
declare -a DLD1_CONTROL_BAMS=()

while IFS=, read -r sample_id client_label cell_line dox_status cohort_role primary_role bam_path bam_gib; do
  [[ "$sample_id" == "sample_id" ]] && continue
  [[ -s "$bam_path" && -s "$bam_path.bai" ]] || { echo "BAM/BAI missing for $sample_id" >&2; exit 2; }
  raw_bam="$RAW/${sample_id}.roi_all_records.bam"
  caller_bam="$CALLER/${sample_id}.roi_gridss_compatible.bam"
  audit_tsv="$AUDIT/${sample_id}.filter_impact.tsv"
  if [[ ! -s "$caller_bam" || ! -s "$caller_bam.bai" || ! -s "$audit_tsv" ]]; then
    [[ ! -e "$raw_bam" ]] || { echo "Incomplete historical raw BAM exists: $raw_bam" >&2; exit 2; }
    "$EXTRACT" --targetbed "$BED" --targetmargin 0 --threads 2 \
      --workingdir "$VALID/gridss_extract_work/$sample_id" --output "$raw_bam" "$bam_path"
    "$PY" "$AUDITOR" --sample-id "$sample_id" --raw-bam "$raw_bam" \
      --caller-bam "$caller_bam" --audit-tsv "$audit_tsv" --threads 2
  fi
  GRIDSS_BAMS+=("$caller_bam")
  case "$primary_role" in
    DLD1_case) DLD1_CASE_BAMS+=("$bam_path") ;;
    DLD1_baseline|DLD1_control) DLD1_CONTROL_BAMS+=("$bam_path") ;;
  esac
done < "$SAMPLES"

[[ ${#GRIDSS_BAMS[@]} -eq 18 ]] || { echo "Expected 18 GRIDSS caller BAMs" >&2; exit 2; }
[[ ${#DLD1_CASE_BAMS[@]} -eq 2 && ${#DLD1_CONTROL_BAMS[@]} -eq 3 ]] || {
  echo "DLD1 case/control design assertion failed" >&2; exit 2;
}

cat "$AUDIT"/*.filter_impact.tsv | awk 'NR == 1 || FNR > 1' > "$AUDIT/gridss_filter_impact_all_samples.tsv"

gridss_vcf="$VALID/gridss_joint/gridss_joint_18samples.vcf.gz"
gridss_assembly="$VALID/gridss_joint/gridss_joint_18samples.assembly.bam"
if [[ ! -s "$gridss_vcf" ]]; then
  "$GRIDSS" -r "$GRIDSS_REF" -o "$gridss_vcf" -a "$gridss_assembly" \
    -t 20 --jvmheap 48g --otherjvmheap 4g \
    -w "$VALID/gridss_joint/work" "${GRIDSS_BAMS[@]}"
fi
[[ -s "$gridss_vcf" && -s "$gridss_assembly" ]] || { echo "GRIDSS output missing" >&2; exit 2; }

svaba_prefix="$VALID/svaba_dld1/svaba_dld1_targeted"
if [[ ! -s "${svaba_prefix}.bps.txt.gz" ]]; then
  (
    cd "$VALID/svaba_dld1"
    "$SVABA" run -t "${DLD1_CASE_BAMS[@]}" -n "${DLD1_CONTROL_BAMS[@]}" \
      -G "$REF" -a "$(basename "$svaba_prefix")" -p 8 -k "$BED" \
      --read-tracking --all-contigs -z
  )
fi
[[ -s "${svaba_prefix}.bps.txt.gz" ]] || { echo "SvABA breakpoint output missing" >&2; exit 2; }

touch "$VALID/.complete"
echo "EXTERNAL_SV_VALIDATION_COMPLETE $(date -Is)" | tee -a "$RUNLOG"
