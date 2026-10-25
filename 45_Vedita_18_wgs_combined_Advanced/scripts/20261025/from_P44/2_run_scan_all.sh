#!/usr/bin/env bash
# Runs Stage 2 for all 18 samples, 6 at a time (6 x 4 threads = 24 threads; I/O bound).
# Resumable: finished samples carry stage2_scan/<S>/.done and are skipped.
source "$(dirname "$0")/0_common.sh"
samples_all | xargs -P 6 -I{} bash "$PROJ/scripts/2_genome_scan_rescue.sh" {} > "$LOG/stage2_scan.log" 2>&1
n=$(ls "$OUT"/stage2_scan/*/.done 2>/dev/null | wc -l)
echo "stage2 finished: $n/18 samples done" >> "$LOG/stage2_scan.log"
[[ $n -eq 18 ]]
