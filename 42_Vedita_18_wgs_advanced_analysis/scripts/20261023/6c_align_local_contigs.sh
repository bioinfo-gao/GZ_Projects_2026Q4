#!/usr/bin/env bash
# Combine successful P42 local-assembly contigs and align every placement to GRCh38.
# `-a -Y` intentionally retains alternative/supplementary placements for repeat-aware review.
set -euo pipefail

ASSEMBLY_DIR=${1:?"Usage: $0 <6b_assembly_dir> <6c_output_dir>"}
OUT_DIR=${2:?"Usage: $0 <6b_assembly_dir> <6c_output_dir>"}
BWA=/Work_bio/gao/configs/.conda/envs/regular_bioinfo/bin/bwa-mem2
PY=/Work_bio/gao/configs/.conda/envs/regular_bioinfo/bin/python
INDEX=/Work_bio/references/Homo_sapiens/GRCh38/GATK.GRCh38/BWAmem2Index/Homo_sapiens_assembly38.fasta

[[ -f "$ASSEMBLY_DIR/.complete" && -s "$ASSEMBLY_DIR/assembly_run_status.tsv" ]] || {
  echo "Completed 6b assembly is required: $ASSEMBLY_DIR" >&2; exit 2;
}
[[ -x "$BWA" && -s "${INDEX}.0123" ]] || { echo "Missing bwa-mem2 or index: $INDEX" >&2; exit 2; }
[[ ! -e "$OUT_DIR" ]] || { echo "Refusing to overwrite contig-alignment output: $OUT_DIR" >&2; exit 2; }
mkdir -p "$OUT_DIR"

combined="$OUT_DIR/all_local_assembly_contigs.fasta"
: > "$combined"
while IFS=$'\t' read -r locus_id sample_id paired unpaired recovered status detail; do
  [[ "$locus_id" == locus_id || "$status" != completed ]] && continue
  fasta="$ASSEMBLY_DIR/loci/$locus_id/megahit/final.contigs.fa"
  [[ -s "$fasta" ]] || { echo "Completed locus lacks contig FASTA: $locus_id" >&2; exit 3; }
  awk -v prefix="${locus_id}" '/^>/ {sub(/^>/, ">" prefix "|"); print; next} {print}' "$fasta" >> "$combined"
done < "$ASSEMBLY_DIR/assembly_run_status.tsv"
[[ -s "$combined" ]] || { echo "No contigs available for placement alignment" >&2; exit 3; }

"$BWA" mem -a -Y -t 20 "$INDEX" "$combined" > "$OUT_DIR/all_local_assembly_contigs.all_placements.sam" 2> "$OUT_DIR/bwa_mem2.log"
"$PY" "$PWD/scripts/20261023/6d_summarize_contig_placements.py" \
  --contigs "$combined" --sam "$OUT_DIR/all_local_assembly_contigs.all_placements.sam" --out-dir "$OUT_DIR"
printf 'completed\n' > "$OUT_DIR/.complete"
