#!/usr/bin/env bash
# Indexed local depth at every base in the 121.5-kb discovery window.
# Both tracks use primary, non-duplicate, non-QC-fail, non-supplementary records
# and suppress mate-overlap double counting (-s); only MAPQ differs.
set -euo pipefail
source "$(dirname "$0")/0_common.sh"
STAGE="$OUT/03_coverage/discovery_per_base"
mkdir -p "$STAGE"
failed=0; pids=()

run_one() {
  local sid="$1" bam="$2" sdir="$STAGE/$sid"
  mkdir -p "$sdir"
  if [[ -f "$sdir/.complete" ]]; then
    echo "[$(date -Is)] $sid already complete; preserving existing result"
    return 0
  fi
  "$SAMTOOLS" depth -aa -s -r "$DISCOVERY_REGION" -G 2048 -Q 20 "$bam" \
    | gzip -c > "$sdir/standard_primary_mapq20.depth.tsv.gz"
  "$SAMTOOLS" depth -aa -s -r "$DISCOVERY_REGION" -G 2048 -Q 0 "$bam" \
    | gzip -c > "$sdir/inclusive_mapq0.depth.tsv.gz"
  test -s "$sdir/standard_primary_mapq20.depth.tsv.gz"
  test -s "$sdir/inclusive_mapq0.depth.tsv.gz"
  printf 'completed %s\n' "$(date -Is)" > "$sdir/.complete"
}
wait_batch() { local pid; for pid in "${pids[@]}"; do wait "$pid" || failed=1; done; pids=(); }
while IFS=, read -r sid _ _ _ _ _ bam _; do
  [[ "$sid" == sample_id ]] && continue
  run_one "$sid" "$bam" >"$LOGDIR/03b_discovery_depth_${sid}.log" 2>&1 & pids+=("$!")
  [[ ${#pids[@]} -eq 4 ]] && wait_batch
done < "$SAMPLES"
wait_batch
[[ "$failed" -eq 0 ]] || { echo 'Local depth job failed; inspect per-sample logs.' >&2; exit 1; }
echo "[$(date -Is)] discovery per-base depth complete: $STAGE"
