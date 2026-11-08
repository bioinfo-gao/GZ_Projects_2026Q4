#!/usr/bin/env bash
# P51 — resident auto-chain (tmux p51_watch). Every 10 min:
#   1. refresh sample table; are all 24 samples paired on disk?
#   2. partial phase finished (PIPELINE_COMPLETE_partial.marker) and p51_rnaseq gone?
#      -> launch `2_nextflow.sh full` (-resume reuses all partial per-sample tasks)
#   3. full phase finished -> run 7_orchestrate_downstream.sh (DE -> enrichment -> report
#      -> self-audit -> delivery folder), once.
# Never kills a running pipeline. Stops itself on any FAILED marker.
# Usage: bash scripts/3_wait_and_launch.sh     (self-relaunches into tmux p51_watch)
set -uo pipefail
SCRIPT="$(readlink -f "${BASH_SOURCE[0]}")"
SCRIPT_DIR="$(dirname "$SCRIPT")"
PROJ_DIR="$(dirname "$SCRIPT_DIR")"
L="${PROJ_DIR}/logs"; mkdir -p "$L"
WLOG="${L}/watch_$(date '+%Y%m%d').log"

if [ -z "${TMUX:-}" ]; then
    tmux has-session -t p51_watch 2>/dev/null && { echo "p51_watch already running"; exit 1; }
    tmux new-session -d -s p51_watch "bash '${SCRIPT}' 2>&1 | tee -a '${WLOG}'"
    echo "Launched tmux p51_watch; log: ${WLOG}"; exit 0
fi

PY=/Work_bio/gao/configs/.conda/envs/regular_bioinfo/bin/python
ts() { date '+%F %T'; }
while true; do
    if ls "$L"/PIPELINE_FAILED_*.marker >/dev/null 2>&1; then
        echo "[$(ts)] FAILED marker present -> watcher stops (needs attention)"; exit 1
    fi
    [ -f "$L/DOWNSTREAM_COMPLETE.marker" ] && { echo "[$(ts)] downstream complete -> exit"; exit 0; }
    npair=$(cd "$SCRIPT_DIR" && "$PY" 0_build_sample_table.py 2>/dev/null | sed -n 's/.*paired on disk: \([0-9]*\)\/24.*/\1/p')
    running=$(tmux has-session -t p51_rnaseq 2>/dev/null && echo yes || echo no)
    echo "[$(ts)] paired=${npair:-?}/24 pipeline_tmux=${running}"
    if [ "$running" = "no" ]; then
        if [ -f "$L/PIPELINE_COMPLETE_full.marker" ]; then
            if [ ! -f "$L/DOWNSTREAM_LAUNCHED.marker" ]; then
                if [ -f "$SCRIPT_DIR/7_orchestrate_downstream.sh" ]; then
                    touch "$L/DOWNSTREAM_LAUNCHED.marker"
                    echo "[$(ts)] full pipeline done -> launching downstream"
                    bash "$SCRIPT_DIR/7_orchestrate_downstream.sh"
                else
                    echo "[$(ts)] full done but 7_orchestrate_downstream.sh missing -> waiting"
                fi
            fi
        elif [ -f "$L/PIPELINE_COMPLETE_partial.marker" ] && [ "${npair:-0}" = "24" ]; then
            echo "[$(ts)] 24/24 on disk + partial done -> launching FULL"
            bash "$SCRIPT_DIR/2_nextflow.sh" full
        fi
    fi
    sleep 600
done
