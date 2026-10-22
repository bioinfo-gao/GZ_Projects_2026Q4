#!/usr/bin/env bash
# P42 formal read-rescue workflow.  It writes only a new dated analysis directory.
set -euo pipefail
SCRIPT="$(readlink -f "$0")"
source "$(dirname "$SCRIPT")/0_common.sh"
LOG="$PROJ/logs/P42_read_rescue_${RUN_DATE}.log"
if [[ -z "${TMUX:-}" ]]; then
  mkdir -p "$PROJ/logs"
  tmux new-session -d -s p42_read_rescue "bash '$SCRIPT' 2>&1 | tee '$LOG'"
  echo "Started tmux session p42_read_rescue; log: $LOG"
  exit 0
fi
need "$REF"; need "$BWA_REF.0123"; need "$SAMPLES"; need "$BWA"
"$PY" "$SCRIPTS/1_prepare_candidates.py" "$SAMPLES" "$P41_EVIDENCE" "$RUN/candidates"
"$PY" "$SCRIPTS/2_extract_support_reads.py" "$SAMPLES" "$RUN/candidates/candidate_support_fragments.tsv" "$WINDOW" "$RUN/reads"
mkdir -p "$RUN/remap"
for fq in "$RUN"/reads/S*.fastq; do
  sample=$(basename "$fq" .fastq)
  "$BWA" mem -a -Y -t 8 "$BWA_REF" "$fq" > "$RUN/remap/$sample.sam"
done
first=1
for sam in "$RUN"/remap/S*.sam; do
  if [[ "$first" -eq 1 ]]; then
    cp "$sam" "$RUN/remap/all_samples.sam"
    first=0
  else
    grep -v '^@' "$sam" >> "$RUN/remap/all_samples.sam"
  fi
done
"$PY" "$SCRIPTS/3_classify_rescue.py" "$RUN/candidates/candidate_junctions.tsv" "$RUN/reads/rescue_read_manifest.tsv" "$RUN/remap/all_samples.sam" "$RUN/results"
[[ $(($(wc -l < "$RUN/results/rescue_candidate_results.tsv") - 1)) -eq 179 ]] || die "candidate output row count is not 179"
echo "P42_READ_RESCUE_DONE candidates=179"
