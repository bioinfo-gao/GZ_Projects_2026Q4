#!/usr/bin/env bash
# Run independent MEGAHIT assemblies for P42's prepared sample/local-locus inputs.
# Default concurrency is intentionally bounded at 4 x 6 threads = 24 CPU threads;
# each locus has its own resumable inputs, log, status, and output directory.
set -euo pipefail

PROJECT=/home/gao/projects_2026Q4/42_Vedita_18_wgs_advanced_analysis
INPUT_DIR=${1:?"Usage: $0 <6a_inputs_dir> <assembly_out_dir> [threads_per_job] [max_parallel_jobs]"}
OUT_DIR=${2:?"Usage: $0 <6a_inputs_dir> <assembly_out_dir> [threads_per_job] [max_parallel_jobs]"}
THREADS_PER_JOB=${3:-6}
MAX_PARALLEL_JOBS=${4:-4}
MEGAHIT=/Work_bio/gao/configs/.conda/envs/regular_bioinfo/bin/megahit

LOCI_TSV="$INPUT_DIR/assembly_loci.tsv"
[[ -s "$LOCI_TSV" && -f "$INPUT_DIR/.complete" ]] || { echo "6a inputs are incomplete: $INPUT_DIR" >&2; exit 2; }
[[ -f "$INPUT_DIR/unmapped_mate_rescue/.complete" ]] || { echo "Unmapped-mate rescue is incomplete: $INPUT_DIR/unmapped_mate_rescue" >&2; exit 2; }
[[ -x "$MEGAHIT" ]] || { echo "Missing MEGAHIT: $MEGAHIT" >&2; exit 2; }
[[ $((THREADS_PER_JOB * MAX_PARALLEL_JOBS)) -le 28 ]] || {
  echo "Thread budget exceeds P42 physical-core limit: ${THREADS_PER_JOB} x ${MAX_PARALLEL_JOBS}" >&2; exit 2;
}
[[ ! -e "$OUT_DIR" ]] || { echo "Refusing to overwrite assembly output: $OUT_DIR" >&2; exit 2; }
mkdir -p "$OUT_DIR/logs" "$OUT_DIR/loci"

run_one() {
  local locus_id=$1 sample_id=$2
  local input_locus="$INPUT_DIR/loci/$locus_id"
  local output_locus="$OUT_DIR/loci/$locus_id"
  local r1="$input_locus/reads_R1.fastq.gz" r2="$input_locus/reads_R2.fastq.gz" single="$input_locus/reads_unpaired.fastq.gz"
  local recovered="$INPUT_DIR/unmapped_mate_rescue/loci/$locus_id/reads_recovered_mates.fastq.gz"
  local paired_count unpaired_count recovered_count single_inputs status=completed detail=success
  paired_count=$(gzip -cd "$r1" | awk 'END {print NR/4}')
  unpaired_count=$(gzip -cd "$single" | awk 'END {print NR/4}')
  recovered_count=0
  single_inputs="$single"
  if [[ -s "$recovered" ]]; then
    recovered_count=$(gzip -cd "$recovered" | awk 'END {print NR/4}')
    single_inputs="$single,$recovered"
  fi
  mkdir -p "$output_locus"
  if (( paired_count < 2 && unpaired_count < 3 )); then
    status=not_run_insufficient_sequences
    detail="fewer_than_2_paired_templates_and_3_unpaired_read_ends"
  else
    if ! "$MEGAHIT" -1 "$r1" -2 "$r2" -r "$single_inputs" \
        -t "$THREADS_PER_JOB" -m 12000000000 \
        --k-list 21,31,41,51,61,71,81,91,101,111,121,131,141 \
        --min-contig-len 100 -o "$output_locus/megahit" \
        > "$OUT_DIR/logs/${locus_id}.megahit.log" 2>&1; then
      status=failed
      detail=megahit_nonzero_exit
    elif [[ ! -s "$output_locus/megahit/final.contigs.fa" ]]; then
      status=failed
      detail=missing_or_empty_final_contigs
    fi
  fi
  printf 'locus_id\tsample_id\tpaired_templates\tunpaired_read_ends\trecovered_unmapped_read_ends\tstatus\tdetail\n%s\t%s\t%s\t%s\t%s\t%s\t%s\n' \
    "$locus_id" "$sample_id" "$paired_count" "$unpaired_count" "$recovered_count" "$status" "$detail" > "$output_locus/status.tsv"
  [[ "$status" != failed ]]
}

pids=() ; names=()
while IFS=$'\t' read -r locus_id sample_id _; do
  [[ "$locus_id" == locus_id ]] && continue
  run_one "$locus_id" "$sample_id" &
  pids+=("$!"); names+=("$locus_id")
  if (( ${#pids[@]} >= MAX_PARALLEL_JOBS )); then
    for i in "${!pids[@]}"; do
      wait "${pids[$i]}" || echo "FAILED ${names[$i]}" >&2
    done
    pids=(); names=()
  fi
done < "$LOCI_TSV"
for i in "${!pids[@]}"; do
  wait "${pids[$i]}" || echo "FAILED ${names[$i]}" >&2
done

status_files=$(find "$OUT_DIR/loci" -mindepth 2 -maxdepth 2 -type f -name status.tsv | sort)
status_count=$(printf '%s\n' "$status_files" | sed '/^$/d' | wc -l)
expected_count=$(($(wc -l < "$LOCI_TSV") - 1))
[[ "$status_count" -eq "$expected_count" ]] || { echo "Expected $expected_count locus status files, found $status_count" >&2; exit 3; }
{
  head -1 "$(printf '%s\n' "$status_files" | head -1)"
  while IFS= read -r status_file; do tail -n +2 "$status_file"; done <<< "$status_files"
} > "$OUT_DIR/assembly_run_status.tsv"
failed_count=$(awk -F'\t' 'NR>1 && $6=="failed" {n++} END {print n+0}' "$OUT_DIR/assembly_run_status.tsv")
not_run_count=$(awk -F'\t' 'NR>1 && $6=="not_run_insufficient_sequences" {n++} END {print n+0}' "$OUT_DIR/assembly_run_status.tsv")
printf '{\n  "input_dir": "%s",\n  "threads_per_job": %s,\n  "max_parallel_jobs": %s,\n  "physical_thread_budget": %s,\n  "expected_loci": %s,\n  "failed_loci": %s,\n  "not_run_insufficient_sequences": %s\n}\n' \
  "$INPUT_DIR" "$THREADS_PER_JOB" "$MAX_PARALLEL_JOBS" "$((THREADS_PER_JOB * MAX_PARALLEL_JOBS))" "$expected_count" "$failed_count" "$not_run_count" > "$OUT_DIR/audit.json"
(( failed_count == 0 )) || exit 4
printf 'completed\n' > "$OUT_DIR/.complete"
