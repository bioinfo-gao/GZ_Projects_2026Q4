#!/usr/bin/env bash
# Standard CN track: primary, non-duplicate, non-QC-fail, non-supplementary, MAPQ >=20.
# Four BAM scans run concurrently in tmux; each mosdepth process uses two threads.
set -euo pipefail
source "$(dirname "$0")/0_common.sh"

STAGE="$OUT/03_coverage/standard_genome100kb"
BEDDIR="$OUT/03_coverage/bins"
mkdir -p "$STAGE" "$BEDDIR"
"$PY" "$(dirname "$0")/2_make_coverage_bins.py" --fai "${REF}.fai" --out-dir "$BEDDIR"

failed=0; pids=()
run_one() {
  local sid="$1" bam="$2" prefix="$STAGE/$sid/$sid.standard_mapq20_primary"
  mkdir -p "$STAGE/$sid"
  if [[ -f "${prefix}.regions.bed.gz" && -f "$STAGE/$sid/.complete" ]]; then
    echo "[$(date -Is)] $sid already complete; preserving existing result"
    return 0
  fi
  "$MOSDEPTH" --threads 2 --no-per-base --by "$BEDDIR/canonical_100kb.bed" \
    --flag 3844 --mapq 20 "$prefix" "$bam"
  test -s "${prefix}.regions.bed.gz"
  printf 'completed %s\n' "$(date -Is)" > "$STAGE/$sid/.complete"
}
wait_batch() { local pid; for pid in "${pids[@]}"; do wait "$pid" || failed=1; done; pids=(); }
while IFS=, read -r sid _ _ _ _ _ bam _; do
  [[ "$sid" == sample_id ]] && continue
  run_one "$sid" "$bam" >"$LOGDIR/03_coverage_${sid}.log" 2>&1 & pids+=("$!")
  [[ ${#pids[@]} -eq 4 ]] && wait_batch
done < "$SAMPLES"
wait_batch
[[ "$failed" -eq 0 ]] || { echo 'Coverage job failed; inspect per-sample logs.' >&2; exit 1; }
echo "[$(date -Is)] standard genome coverage complete: $STAGE"
