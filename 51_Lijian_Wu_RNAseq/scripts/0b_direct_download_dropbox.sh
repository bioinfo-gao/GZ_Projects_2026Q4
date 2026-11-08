#!/usr/bin/env bash
# P51 — fallback when the local Dropbox client stalls (2026-11-08 01:2x: daemon reported
# "Downloading 78 files" but had received ~0 B/s since 22:51; the daemon itself was NOT touched).
# Fetches the missing FASTQs straight from the Dropbox API (single-use temporary links issued by
# the Dropbox connector) into a SEPARATE dir; never writes into the Dropbox-synced folder.
# Usage: bash 0b_direct_download_dropbox.sh <links.tsv>   (cols: name  size  content_hash  url)
# Each file -> <name>.part, renamed only after size == cloud size; content_hash checked by 0c.
set -uo pipefail
LINKS="$1"
OUT=/home/gao/projects_2026Q4/51_Lijian_Wu_RNAseq/raw_fastq_direct_dropbox_20261108
mkdir -p "$OUT"
while IFS=$'\t' read -r name size hash url; do
  (
    if [ -f "$OUT/$name" ] && [ "$(stat -c%s "$OUT/$name")" = "$size" ]; then echo "SKIP $name"; exit 0; fi
    curl -sS --fail -L -o "$OUT/$name.part" "$url"; rc=$?
    got=$(stat -c%s "$OUT/$name.part" 2>/dev/null || echo 0)
    if [ $rc -eq 0 ] && [ "$got" = "$size" ]; then mv "$OUT/$name.part" "$OUT/$name"; echo "$(date '+%F %T') DONE $name"
    else echo "$(date '+%F %T') FAIL $name rc=$rc got=$got want=$size"; fi
  ) &
done < "$LINKS"
wait
echo "$(date '+%F %T') ALL_FINISHED"
