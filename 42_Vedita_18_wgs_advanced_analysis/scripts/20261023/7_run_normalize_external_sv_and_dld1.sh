#!/usr/bin/env bash
# P42 stage 7: normalize completed external caller outputs and DLD1 matrix.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/0_common.sh"

OUTDIR="$OUT/07_external_sv_dld1_normalization_v1"
GRIDSS="$OUT/06_external_sv_validation_v3/gridss_joint/gridss_joint_18samples.vcf.gz"
SVABA_DIR="$OUT/06_external_svaba_dld1_v5_multisample/svaba_dld1"
CLUSTERS="$OUT/04_junction_evidence/clusters/candidate_clusters.tsv"
CONTIGS="$OUT/05_local_assembly/contig_placements_v1/contig_placement_summary.tsv"

[[ ! -e "$OUTDIR" ]] || { echo "Refusing to overwrite existing stage: $OUTDIR" >&2; exit 2; }
for f in "$GRIDSS" "$SVABA_DIR/svaba_dld1_targeted_multisample.svaba.germline.sv.vcf.gz" \
  "$SVABA_DIR/svaba_dld1_targeted_multisample.svaba.somatic.sv.vcf.gz" "$CLUSTERS" "$CONTIGS"; do
  [[ -s "$f" ]] || { echo "Required input missing: $f" >&2; exit 2; }
done

"$PY" "$SCRIPT_DIR/7_normalize_external_sv_and_dld1.py" \
  --gridss-vcf "$GRIDSS" \
  --svaba-germline-vcf "$SVABA_DIR/svaba_dld1_targeted_multisample.svaba.germline.sv.vcf.gz" \
  --svaba-somatic-vcf "$SVABA_DIR/svaba_dld1_targeted_multisample.svaba.somatic.sv.vcf.gz" \
  --cluster-tsv "$CLUSTERS" \
  --contig-summary "$CONTIGS" \
  --out-dir "$OUTDIR"

[[ -s "$OUTDIR/dld1_external_caller_matrix.tsv" && -s "$OUTDIR/audit.json" && -s "$OUTDIR/.complete" ]] || {
  echo "Stage 7 outputs incomplete" >&2; exit 2;
}
