#!/usr/bin/env bash
# P51 — nf-core/rnaseq 3.15.1, mouse GRCm39 / GENCODE M35, solo full-machine profile.
#
# Usage:  bash scripts/2_nextflow.sh partial   # whatever samples are complete on disk now
#         bash scripts/2_nextflow.sh full      # refuses to run unless all 24 samples pass
#
# Dropbox is still syncing (2026-11-07 22:49: 13/24 paired). Phase "partial" processes the
# samples already here; phase "full" re-runs with the complete samplesheet and -resume, so
# every per-sample task already done in "partial" is reused from work/ cache (same inputs,
# same params, same work dir). Downstream DE is NOT run on partial data.
#
# Self-relaunches into tmux session p51_rnaseq; auto-retries once with -resume.
set -uo pipefail

PHASE="${1:-partial}"
SCRIPT="$(readlink -f "${BASH_SOURCE[0]}")"
SCRIPT_DIR="$(dirname "$SCRIPT")"
PROJ_DIR="$(dirname "$SCRIPT_DIR")"
STAMP="$(date '+%Y%m%d_%H%M')"
mkdir -p "${PROJ_DIR}/logs"                      # tee dies silently if its dir is missing
LOG="${PROJ_DIR}/logs/nextflow_${PHASE}_${STAMP}.log"

if [ -z "${TMUX:-}" ]; then
    if tmux has-session -t p51_rnaseq 2>/dev/null; then
        echo "tmux session p51_rnaseq already exists — refusing to start a second run"; exit 1
    fi
    tmux new-session -d -s p51_rnaseq "bash '${SCRIPT}' '${PHASE}' 2>&1 | tee '${LOG}'"
    echo "Launched tmux p51_rnaseq (phase=${PHASE}); log: ${LOG}"
    exit 0
fi

cd "$SCRIPT_DIR" || exit 1
PY=/Work_bio/gao/configs/.conda/envs/regular_bioinfo/bin/python
NEXTFLOW=/Work_bio/gao/configs/.conda/envs/regular_bioinfo/bin/nextflow
REF_DIR=/Work_bio/references/Mus_musculus/GRCm39/mouse_gencode_M35
export NXF_SINGULARITY_CACHEDIR="${NXF_SINGULARITY_CACHEDIR:-/home/gao/.singularity/nf-core}"
export NXF_OPTS="-Xms512m -Xmx4g"

echo "=== [$(date '+%F %T')] phase=${PHASE}: refresh sample table + samplesheet ==="
"$PY" 0_build_sample_table.py || exit 1
if [ "$PHASE" = "full" ]; then
    "$PY" 1_produce_nf-core_Samplesheet.py --require-all || { echo "NOT ALL 24 READY"; exit 1; }
else
    "$PY" 1_produce_nf-core_Samplesheet.py || exit 1
fi
cp nf_core_samplesheet.csv "${PROJ_DIR}/logs/nf_core_samplesheet_${PHASE}_${STAMP}.csv"

run_nextflow() {
    "$NEXTFLOW" run nf-core/rnaseq \
        -r 3.15.1 \
        -profile singularity \
        -c local_optimized.config \
        -resume \
        --input  "${SCRIPT_DIR}/nf_core_samplesheet.csv" \
        --outdir "${PROJ_DIR}/output_results" \
        --fasta      "${REF_DIR}/GRCm39.primary_assembly.genome.fa" \
        --gtf        "${REF_DIR}/gencode.vM35.annotation.gtf" \
        --star_index "${REF_DIR}/star_index" \
        --extra_star_align_args '--twopassMode None --outFilterMultimapNmax 3' \
        --gencode \
        --aligner star_salmon \
        --max_cpus 56 \
        --max_memory '108.GB' \
        -work-dir "${PROJ_DIR}/work"
}

echo "=== [$(date '+%F %T')] nf-core/rnaseq 3.15.1 phase=${PHASE} attempt 1 ==="
if run_nextflow; then
    echo "=== [$(date '+%F %T')] PIPELINE_COMPLETE phase=${PHASE} ==="
    touch "${PROJ_DIR}/logs/PIPELINE_COMPLETE_${PHASE}.marker"; exit 0
fi
echo "=== [$(date '+%F %T')] attempt 1 failed — retry with -resume ==="
if run_nextflow; then
    echo "=== [$(date '+%F %T')] PIPELINE_COMPLETE phase=${PHASE} (retry) ==="
    touch "${PROJ_DIR}/logs/PIPELINE_COMPLETE_${PHASE}.marker"; exit 0
fi
echo "=== [$(date '+%F %T')] PIPELINE_FAILED phase=${PHASE} — human attention needed ==="
touch "${PROJ_DIR}/logs/PIPELINE_FAILED_${PHASE}.marker"
exit 1
