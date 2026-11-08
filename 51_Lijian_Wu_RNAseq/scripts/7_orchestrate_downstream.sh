#!/usr/bin/env bash
# P51 — downstream chain after the FULL nf-core run (called by 3_wait_and_launch.sh, or by hand).
#   4a fill read pairs -> 4 DESeq2/PCA -> 5 enrichment -> READY_FOR_REPORT.marker
# Report writing (6), pre-delivery audit (8) and delivery folder (9) follow once results
# have been read by the analyst — the report's Core Conclusions interpret the biology.
# Logs: logs/downstream_<stamp>_<step>.log ; Rscript via absolute path, `> log 2>&1`
# (tee would mask the Rscript exit code).
set -uo pipefail
S="$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")"; P="$(dirname "$S")"; L="$P/logs"
STAMP="$(date '+%Y%m%d_%H%M')"
PY=/Work_bio/gao/configs/.conda/envs/regular_bioinfo/bin/python
RS=/Work_bio/gao/configs/.conda/envs/DE_R45/bin/Rscript
step() { echo "[$(date '+%F %T')] >>> $1"; }
fail() { echo "[$(date '+%F %T')] !!! $1 FAILED (see $2)"; touch "$L/DOWNSTREAM_FAILED_$1.marker"; exit 1; }
cd "$S" || exit 1
step 4a; "$PY" 4a_fill_read_pairs.py > "$L/downstream_${STAMP}_4a.log" 2>&1 || fail 4a "$L/downstream_${STAMP}_4a.log"
step 4;  "$RS" 4_run_DE_PCA.R        > "$L/downstream_${STAMP}_4.log"  2>&1 || fail 4  "$L/downstream_${STAMP}_4.log"
step 5;  "$RS" 5_run_enrichment.R    > "$L/downstream_${STAMP}_5.log"  2>&1 || fail 5  "$L/downstream_${STAMP}_5.log"
touch "$L/READY_FOR_REPORT.marker"
step "DE + enrichment done -> READY_FOR_REPORT (analyst writes report, audits, builds delivery)"
