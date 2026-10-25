#!/usr/bin/env bash
# P45 analysis chain after import: waits for Stage 0 (.done), then 8x -> 9x -> 10 -> 10c -> 11 -> 14 -> 14b.
# Stops at the first failing stage (set -e). Every stage writes only under analysis/20261025_p45/ (10c: docs/).
source "$(dirname "$0")/0_common.sh"
cd "$SCRIPTS"
until [[ -f $INH/.done ]]; do
  tmux has-session -t p45_stage0_import 2>/dev/null || { echo "Stage 0 ended without .done"; exit 2; }
  sleep 60
done
[[ -f $INH/P42/external_callers/.done ]] || { echo "Stage 0b missing"; exit 1; }
PRIMER_PY=$ASM_ENV/python            # primer3-py lives in the P44 assembly env
for st in 8x_caller_concordance_4callers 9x_structural_coherence 10_build_tables 10c_p42_reconciliation 11_figures \
          14_primer_design 14b_primer_locus_check; do
  echo "$(date '+%F %T') $st start"
  py=$PY; [[ $st == 14_primer_design ]] && py=$PRIMER_PY
  $py -u "$SCRIPTS/$st.py" > "$LOG/$st.log" 2>&1 || { echo "$st FAILED"; tail -20 "$LOG/$st.log"; exit 1; }
  tail -1 "$LOG/$st.log"
done
echo "$(date '+%F %T') P45_CHAIN_DONE"
