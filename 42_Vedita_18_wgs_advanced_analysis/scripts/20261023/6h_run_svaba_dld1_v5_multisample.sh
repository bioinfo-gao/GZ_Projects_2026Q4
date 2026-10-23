#!/usr/bin/env bash
# P42 v5: all five DLD1 BAMs must be passed with repeated SvABA input flags.
# SvABA accepts multiple BAMs, but they must each retain their own -t/-n option.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/0_common.sh"

VALID="${P42_SVABA_STAGE:-$OUT/06_external_svaba_dld1_v5_multisample}"
GRIDSS_REF="$OUT/06_external_sv_validation/gridss_reference_v1/Homo_sapiens_assembly38.fasta"
SVABA=/Work_bio/gao/configs/.conda/envs/p42_sv_validation_20261023/bin/svaba
BED="$VALID/regions/discovery_window_0based.bed"
PREFIX="$VALID/svaba_dld1/svaba_dld1_targeted_multisample"
MANIFEST="$VALID/svaba_dld1/input_manifest.tsv"

mkdir -p "$VALID/regions" "$VALID/svaba_dld1"
[[ -s "$GRIDSS_REF" && -s "$GRIDSS_REF.bwt" && -s "$GRIDSS_REF.sa" && -s "$GRIDSS_REF.fai" ]] || {
  echo "Project-private reference/BWA index incomplete: $GRIDSS_REF" >&2; exit 2;
}
[[ -x "$SVABA" ]] || { echo "SvABA executable missing: $SVABA" >&2; exit 2; }
printf 'chr11\t114553999\t114675500\tP42_discovery_window\n' > "$BED"

declare -a CASE_BAMS=() CONTROL_BAMS=() CASE_IDS=() CONTROL_IDS=()
while IFS=, read -r sample_id client_label cell_line dox_status cohort_role primary_role bam_path bam_gib; do
  [[ "$sample_id" == "sample_id" ]] && continue
  case "$primary_role" in
    DLD1_case) CASE_BAMS+=("$bam_path"); CASE_IDS+=("$sample_id") ;;
    DLD1_baseline|DLD1_control) CONTROL_BAMS+=("$bam_path"); CONTROL_IDS+=("$sample_id") ;;
    *) continue ;;
  esac
  [[ -s "$bam_path" && -s "$bam_path.bai" ]] || {
    echo "DLD1 BAM/BAI missing for $sample_id: $bam_path" >&2; exit 2;
  }
done < "$SAMPLES"

[[ ${#CASE_BAMS[@]} -eq 2 && ${#CONTROL_BAMS[@]} -eq 3 ]] || {
  echo "DLD1 design assertion failed: expected two cases and three controls" >&2; exit 2;
}
{
  printf 'sample_id\trole\tbam_path\n'
  for i in "${!CASE_BAMS[@]}"; do printf '%s\tcase\t%s\n' "${CASE_IDS[$i]}" "${CASE_BAMS[$i]}"; done
  for i in "${!CONTROL_BAMS[@]}"; do printf '%s\tcontrol\t%s\n' "${CONTROL_IDS[$i]}" "${CONTROL_BAMS[$i]}"; done
} > "$MANIFEST"

if [[ ! -s "${PREFIX}.bps.txt.gz" ]]; then
  declare -a CMD=(run)
  for bam in "${CASE_BAMS[@]}"; do CMD+=(-t "$bam"); done
  for bam in "${CONTROL_BAMS[@]}"; do CMD+=(-n "$bam"); done
  (
    cd "$VALID/svaba_dld1"
    "$SVABA" "${CMD[@]}" -G "$GRIDSS_REF" -a "$(basename "$PREFIX")" -p 8 -k "$BED" \
      --read-tracking --all-contigs -z
  )
fi

[[ -s "${PREFIX}.bps.txt.gz" && -s "${PREFIX}.log" ]] || {
  echo "SvABA outputs missing for multi-sample DLD1 run" >&2; exit 2;
}
[[ $(grep -c '^BAM PARAMS FOR:' "${PREFIX}.log") -eq 5 ]] || {
  echo "SvABA input audit failed: expected exactly 5 parsed BAMs" >&2; exit 2;
}
for bam in "${CASE_BAMS[@]}" "${CONTROL_BAMS[@]}"; do
  grep -F -- "$bam" "${PREFIX}.log" >/dev/null || {
    echo "SvABA input audit failed: missing BAM in tool log: $bam" >&2; exit 2;
  }
done
touch "$VALID/.complete"
echo "SVABA_DLD1_V5_MULTISAMPLE_COMPLETE $(date -Is)"
