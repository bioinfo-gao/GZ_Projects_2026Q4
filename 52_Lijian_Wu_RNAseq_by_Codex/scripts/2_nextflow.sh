#!/usr/bin/env bash
set -euo pipefail
S="$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")"; P="$(dirname "$S")"
if [ -z "${TMUX:-}" ]; then
 tmux new-session -d -s p52_full "bash '$P/scripts/0_run_all.sh' >> '$P/logs/master.log' 2>&1"
 exit 0
fi
cd "$S"
export PATH=/Work_bio/gao/configs/.conda/envs/regular_bioinfo/bin:$PATH
export NXF_SINGULARITY_CACHEDIR=/home/gao/.singularity/nf-core
export NXF_OPTS='-Xms512m -Xmx4g'
run_nextflow() {
 nextflow -log "$P/.nextflow.log" run nf-core/rnaseq -r 3.15.1 -profile singularity -c local_optimized.config "$@" \
 --input "$S/nf_core_samplesheet.csv" --outdir "$P/output_results" \
 --fasta /Work_bio/references/Mus_musculus/GRCm39/mouse_gencode_M35/GRCm39.primary_assembly.genome.fa \
 --gtf /Work_bio/references/Mus_musculus/GRCm39/mouse_gencode_M35/gencode.vM35.annotation.gtf \
 --star_index /Work_bio/references/Mus_musculus/GRCm39/mouse_gencode_M35/star_index \
 --extra_star_align_args '--twopassMode None --outFilterMultimapNmax 3' \
 --gencode --aligner star_salmon --max_cpus 56 --max_memory '108.GB' -work-dir "$P/work"
}
if run_nextflow; then touch "$P/logs/PIPELINE_COMPLETE.marker"; exit 0; fi
run_nextflow -resume
touch "$P/logs/PIPELINE_COMPLETE.marker"
