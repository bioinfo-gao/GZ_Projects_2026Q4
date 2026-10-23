#!/usr/bin/env bash
set -euo pipefail
source "$(dirname "$0")/0_common.sh"
STAGE="$OUT/04_junction_evidence/raw_evidence"
mkdir -p "$STAGE"
failed=0; pids=()
run_one() {
  local sid="$1" bam="$2"
  [[ -f "$STAGE/$sid/.complete" ]] && { echo "[$(date -Is)] $sid already complete"; return; }
  "$PY" -u "$(dirname "$0")/5a_extract_junction_evidence.py" --sample-id "$sid" --bam "$bam" \
    --out-dir "$STAGE/$sid" --chrom "$TARGET_CHROM" --start "$DISCOVERY_START" --end "$DISCOVERY_END"
}
wait_batch() { local pid; for pid in "${pids[@]}"; do wait "$pid" || failed=1; done; pids=(); }
while IFS=, read -r sid _ _ _ _ _ bam _; do
  [[ "$sid" == sample_id ]] && continue
  run_one "$sid" "$bam" >"$LOGDIR/05a_junction_evidence_${sid}.log" 2>&1 & pids+=("$!")
  [[ ${#pids[@]} -eq 4 ]] && wait_batch
done < "$SAMPLES"
wait_batch
[[ "$failed" -eq 0 ]] || { echo 'Evidence extraction failure; inspect logs.' >&2; exit 1; }
echo "[$(date -Is)] raw junction evidence complete: $STAGE"
