#!/usr/bin/env bash
# P42 v4: resume only the SvABA DLD1 validation that v3 could not start.
# v3 GRIDSS outputs and its failed SvABA log remain immutable historical records.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/0_common.sh"

VALID="${P42_SVABA_STAGE:-$OUT/06_external_svaba_dld1_v4}"
GRIDSS_REF="$OUT/06_external_sv_validation/gridss_reference_v1/Homo_sapiens_assembly38.fasta"
SVABA=/Work_bio/gao/configs/.conda/envs/p42_sv_validation_20261023/bin/svaba
BED="$VALID/regions/discovery_window_0based.bed"
PREFIX="$VALID/svaba_dld1/svaba_dld1_targeted"

mkdir -p "$VALID/regions" "$VALID/svaba_dld1"
[[ -s "$GRIDSS_REF" && -s "$GRIDSS_REF.bwt" && -s "$GRIDSS_REF.sa" && -s "$GRIDSS_REF.fai" ]] || {
  echo "Project-private reference/BWA index incomplete: $GRIDSS_REF" >&2; exit 2;
}
[[ -x "$SVABA" ]] || { echo "SvABA executable missing: $SVABA" >&2; exit 2; }

# BED is 0-based, half-open; the source analysis window is 1-based chr11:114554000-114675500.
printf 'chr11\t114553999\t114675500\tP42_discovery_window\n' > "$BED"

declare -a CASE_BAMS=()
declare -a CONTROL_BAMS=()
while IFS=, read -r sample_id client_label cell_line dox_status cohort_role primary_role bam_path bam_gib; do
  [[ "$sample_id" == "sample_id" ]] && continue
  case "$primary_role" in
    DLD1_case) CASE_BAMS+=("$bam_path") ;;
    DLD1_baseline|DLD1_control) CONTROL_BAMS+=("$bam_path") ;;
    *) continue ;;
  esac
  [[ -s "$bam_path" && -s "$bam_path.bai" ]] || {
    echo "DLD1 BAM/BAI missing for $sample_id: $bam_path" >&2; exit 2;
  }
done < "$SAMPLES"

[[ ${#CASE_BAMS[@]} -eq 2 && ${#CONTROL_BAMS[@]} -eq 3 ]] || {
  echo "DLD1 design assertion failed: expected two cases and three controls" >&2; exit 2;
}

if [[ ! -s "${PREFIX}.bps.txt.gz" ]]; then
  (
    cd "$VALID/svaba_dld1"
    "$SVABA" run -t "${CASE_BAMS[@]}" -n "${CONTROL_BAMS[@]}" \
      -G "$GRIDSS_REF" -a "$(basename "$PREFIX")" -p 8 -k "$BED" \
      --read-tracking --all-contigs -z
  )
fi

[[ -s "${PREFIX}.bps.txt.gz" ]] || {
  echo "SvABA breakpoint output missing: ${PREFIX}.bps.txt.gz" >&2; exit 2;
}
touch "$VALID/.complete"
echo "SVABA_DLD1_V4_COMPLETE $(date -Is)"
