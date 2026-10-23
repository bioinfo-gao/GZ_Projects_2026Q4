#!/usr/bin/env bash
# Run four independent sample inventories at a time inside tmux.
set -euo pipefail
source "$(dirname "$0")/0_common.sh"

STAGE="$OUT/01_full_window_inventory"
mkdir -p "$STAGE"
failed=0
pids=()

run_one() {
  local sid="$1"
  if [[ -f "$STAGE/$sid/.complete" ]]; then
    echo "[$(date -Is)] $sid already complete; preserving existing result"
    return 0
  fi
  "$PY" -u "$(dirname "$0")/1_inventory_full_window.py" \
    --sample-table "$SAMPLES" --sample-id "$sid" \
    --chrom "$TARGET_CHROM" --start "$DISCOVERY_START" --end "$DISCOVERY_END" \
    --out-root "$STAGE"
}

wait_batch() {
  local pid
  for pid in "${pids[@]}"; do
    wait "$pid" || failed=1
  done
  pids=()
}

while IFS=, read -r sid _; do
  [[ "$sid" == "sample_id" ]] && continue
  run_one "$sid" >"$LOGDIR/01_inventory_${sid}.log" 2>&1 &
  pids+=("$!")
  if [[ ${#pids[@]} -eq 4 ]]; then
    wait_batch
  fi
done < "$SAMPLES"
wait_batch

if [[ "$failed" -ne 0 ]]; then
  echo "One or more inventory jobs failed; inspect $LOGDIR/01_inventory_*.log" >&2
  exit 1
fi

"$PY" - "$STAGE" <<'PY'
import json, sys
from pathlib import Path
stage = Path(sys.argv[1])
audits = sorted(stage.glob('*/audit.json'))
if len(audits) != 18:
    raise SystemExit(f'Expected 18 audit.json files, found {len(audits)}')
keys = ['records_total', 'records_primary', 'records_secondary',
        'records_supplementary', 'records_mapq0', 'records_with_seq_star',
        'molecules_total', 'unmapped_mate_rescue_attempts',
        'unmapped_mate_rescue_retrieved', 'unmapped_mate_rescue_failures']
out = stage / 'inventory_audit_summary.tsv'
with out.open('w') as h:
    h.write('sample_id\t' + '\t'.join(keys) + '\n')
    for path in audits:
        d = json.loads(path.read_text())
        c = d['record_counts']
        h.write(d['sample_id'] + '\t' + '\t'.join(str(c.get(k, '')) for k in keys) + '\n')
PY

echo "[$(date -Is)] inventory completed: $STAGE/inventory_audit_summary.tsv"
