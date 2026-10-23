#!/usr/bin/env bash
# Persistent P42 stage chain: wait for the active unmapped-mate rescue, then
# assemble 59 local loci and run all-placement contig alignment.  Nonzero exit
# at any stage stops the chain; it never silently skips a failed prerequisite.
set -euo pipefail

PROJECT=/home/gao/projects_2026Q4/42_Vedita_18_wgs_advanced_analysis
INPUTS="$PROJECT/analysis/20261023_fullcohort_repeat_resolve/05_local_assembly/inputs_v1"
ASSEMBLY_OUT="$PROJECT/analysis/20261023_fullcohort_repeat_resolve/05_local_assembly/assemblies_v1"
PLACEMENT_OUT="$PROJECT/analysis/20261023_fullcohort_repeat_resolve/05_local_assembly/contig_placements_v1"
RESCUE_SESSION=p42_unmapped_mate_rescue_20261023

while [[ ! -f "$INPUTS/unmapped_mate_rescue/.complete" ]]; do
  if ! tmux has-session -t "$RESCUE_SESSION" 2>/dev/null; then
    echo "ERROR: unmapped-mate rescue session ended without its completion marker" >&2
    exit 10
  fi
  echo "WAITING: unmapped-mate rescue still active at $(date --iso-8601=seconds)"
  sleep 30
done

echo "START: local assembly at $(date --iso-8601=seconds)"
bash "$PROJECT/scripts/20261023/6b_run_local_assembly.sh" "$INPUTS" "$ASSEMBLY_OUT" 6 4
echo "START: all-placement contig alignment at $(date --iso-8601=seconds)"
bash "$PROJECT/scripts/20261023/6c_align_local_contigs.sh" "$ASSEMBLY_OUT" "$PLACEMENT_OUT"
echo "COMPLETE: P42 local assembly and contig-placement chain at $(date --iso-8601=seconds)"
