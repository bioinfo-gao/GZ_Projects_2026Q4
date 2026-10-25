#!/usr/bin/env bash
# Stage 3b — re-align ONLY the rescued reads (full sequences from primary records) to the same GRCh38 analysis
# set with `bwa-mem2 mem -a -Y`, so every competing placement is written as an explicit secondary record
# (the original BAMs carry them only as XA tags). This is a local re-alignment of the rescue set, not a
# re-alignment of the sequencing run. It tests mapping reproducibility, not biology (plan section 3).
# Usage: 3b_realign_rescue_reads.sh            (all samples whose Stage 3 is done; resumable)
source "$(dirname "$0")/0_common.sh"
for S in $(samples_all); do
  D=$OUT/stage3_rescue/$S
  [[ -f $D/.done ]] || { echo "[$S] stage3 not done — skip"; continue; }
  [[ -f $D/.realign_done ]] && continue
  echo "[$S] $(date '+%T') bwa-mem2 -a"
  $ENV/bwa-mem2 mem -a -Y -t 24 -R "@RG\tID:${S}\tSM:${S}" "$BWA_IDX" \
      "$D/rescue_R1.fq.gz" "$D/rescue_R2.fq.gz" 2> "$D/realign.bwa.log" \
    | $SAMTOOLS sort -@ 4 -o "$D/realigned_a.bam" -
  $SAMTOOLS index "$D/realigned_a.bam"
  touch "$D/.realign_done"
done
