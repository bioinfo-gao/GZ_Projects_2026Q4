#!/usr/bin/env bash
# Chain runner: Stage 3 (6 samples in parallel, each starts as soon as its Stage 2 is done) -> 3b -> 4 -> 5 -> 6 -> 7 -> 8.
# Stops at the first failing stage (set -e); each stage is resumable via its own .done markers / outputs.
source "$(dirname "$0")/0_common.sh"
cd "$PROJ"
stage3_one() {
  local S=$1
  until [[ -f $OUT/stage2_scan/$S/.done ]]; do sleep 30; done
  [[ -f $OUT/stage3_rescue/$S/.done ]] || $PY "$PROJ/scripts/3_collect_rescue_reads.py" "$S" > "$LOG/stage3_$S.log" 2>&1
}
export -f stage3_one; export OUT PROJ PY LOG
echo "$(date '+%F %T') stage3 start"
samples_all | xargs -P 6 -I{} bash -c 'set -e; stage3_one {}'
n=$(ls "$OUT"/stage3_rescue/*/.done | wc -l); [[ $n -eq 18 ]] || { echo "stage3 incomplete: $n/18"; exit 1; }
echo "$(date '+%F %T') stage3b start"; bash "$PROJ/scripts/3b_realign_rescue_reads.sh" > "$LOG/stage3b.log" 2>&1
n=$(ls "$OUT"/stage3_rescue/*/.realign_done | wc -l); [[ $n -eq 18 ]] || { echo "stage3b incomplete: $n/18"; exit 1; }
for st in 4_evidence_clusters 5_local_assembly 6_junction_genotyping 7_breakpoint_annotation 8_caller_concordance; do
  echo "$(date '+%F %T') $st start"
  $PY "$PROJ/scripts/$st.py" > "$LOG/$st.log" 2>&1
  tail -1 "$LOG/$st.log"
done
echo "$(date '+%F %T') CHAIN_DONE"
