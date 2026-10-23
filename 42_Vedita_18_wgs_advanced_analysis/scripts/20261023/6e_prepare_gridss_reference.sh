#!/usr/bin/env bash
# Prepare a project-private traditional-BWA index for GRIDSS without changing
# the shared GRCh38 reference directory. The FASTA dictionary must match a
# representative source BAM before indexing.
set -euo pipefail

REFERENCE=${1:?"Usage: $0 <reference_fasta> <representative_bam> <output_dir>"}
REPRESENTATIVE_BAM=${2:?"Usage: $0 <reference_fasta> <representative_bam> <output_dir>"}
OUT_DIR=${3:?"Usage: $0 <reference_fasta> <representative_bam> <output_dir>"}
SAMTOOLS=/Work_bio/gao/configs/.conda/envs/regular_bioinfo/bin/samtools
BWA=/Work_bio/gao/configs/.conda/envs/p42_sv_validation_20261023/bin/bwa

[[ -s "$REFERENCE" && -s "${REFERENCE}.fai" ]] || { echo "Missing reference FASTA/index: $REFERENCE" >&2; exit 2; }
[[ -s "$REPRESENTATIVE_BAM" && -f "${REPRESENTATIVE_BAM}.bai" ]] || { echo "Missing representative BAM/index: $REPRESENTATIVE_BAM" >&2; exit 2; }
[[ -x "$SAMTOOLS" && -x "$BWA" ]] || { echo "Missing samtools or validation-environment bwa" >&2; exit 2; }
[[ ! -e "$OUT_DIR" ]] || { echo "Refusing to overwrite existing output directory: $OUT_DIR" >&2; exit 2; }

mkdir -p "$OUT_DIR"
local_reference="$OUT_DIR/$(basename "$REFERENCE")"
ln -s "$REFERENCE" "$local_reference"

"$SAMTOOLS" view -H "$REPRESENTATIVE_BAM" | awk -F'\t' '
  $1=="@SQ" {
    sn=""; ln="";
    for (i=2; i<=NF; i++) {
      if ($i ~ /^SN:/) sn=substr($i,4);
      if ($i ~ /^LN:/) ln=substr($i,4);
    }
    print sn "\t" ln;
  }
' > "$OUT_DIR/bam_sequence_dictionary.tsv"
cut -f1,2 "${REFERENCE}.fai" > "$OUT_DIR/reference_sequence_dictionary.tsv"
diff -u "$OUT_DIR/bam_sequence_dictionary.tsv" "$OUT_DIR/reference_sequence_dictionary.tsv" > "$OUT_DIR/bam_reference_dictionary.diff" || {
  echo "BAM/reference sequence dictionaries differ; GRIDSS must not be run" >&2
  exit 3
}
# Keep the zero-byte diff as a machine-readable audit artifact: an empty file
# means the BAM and reference dictionaries matched exactly at preparation time.

"$BWA" index "$local_reference" > "$OUT_DIR/bwa_index.log" 2>&1
for suffix in amb ann bwt pac sa; do
  [[ -s "${local_reference}.${suffix}" ]] || { echo "Missing BWA index component: ${local_reference}.${suffix}" >&2; exit 4; }
done

{
  printf 'field\tvalue\n'
  printf 'source_reference\t%s\n' "$REFERENCE"
  printf 'project_private_reference_link\t%s\n' "$local_reference"
  printf 'representative_bam\t%s\n' "$REPRESENTATIVE_BAM"
  printf 'dictionary_match\tyes\n'
  printf 'bwa\t%s\n' "$($BWA 2>&1 | awk 'NR==2 {print; exit}')"
} > "$OUT_DIR/reference_audit.tsv"
printf 'completed\n' > "$OUT_DIR/.complete"
