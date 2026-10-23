#!/usr/bin/env bash
# P42 stage 9: raw-evidence re-genotyping of every fixed external caller BND.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/0_common.sh"

STAGE7="$OUT/07_external_sv_dld1_normalization_v1"
OUTDIR="$OUT/09_fixed_candidate_regenotyping_v1"
RAW="$OUT/04_junction_evidence/raw_evidence"
COVERAGE="$OUT/03_coverage/discovery_1kb_summary.tsv"
CONTIGS="$OUT/05_local_assembly/contig_placements_v1/contig_placement_summary.tsv"

[[ -s "$STAGE7/.complete" ]] || { echo "Stage 7 is incomplete: $STAGE7" >&2; exit 2; }
[[ ! -e "$OUTDIR" ]] || { echo "Refusing to overwrite existing stage: $OUTDIR" >&2; exit 2; }
for path in "$STAGE7/normalized_external_sv_events.tsv" "$COVERAGE" "$CONTIGS"; do
  [[ -s "$path" ]] || { echo "Required input missing: $path" >&2; exit 2; }
done
for sample in S01_DLD1_pop S02_DLD1_c1_ctrl S03_DLD1_c2_ctrl S04_DLD1_c3_dox S05_DLD1_c4_dox; do
  [[ -s "$RAW/$sample/junction_evidence_records.tsv.gz" ]] || { echo "Raw evidence missing: $sample" >&2; exit 2; }
done

"$PY" "$SCRIPT_DIR/9_regenotype_fixed_candidates.py" \
  --stage7-events "$STAGE7/normalized_external_sv_events.tsv" \
  --raw-evidence-dir "$RAW" \
  --coverage "$COVERAGE" \
  --contigs "$CONTIGS" \
  --out-dir "$OUTDIR"

[[ -s "$OUTDIR/.complete" && -s "$OUTDIR/fixed_candidate_dld1_matrix.tsv" \
  && -s "$OUTDIR/fixed_candidate_summary.tsv" && -s "$OUTDIR/raw_evidence_match_audit.tsv.gz" ]] || {
  echo "Stage 9 outputs incomplete" >&2; exit 2;
}
