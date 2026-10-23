#!/usr/bin/env python3
"""Extract the original full sequence for every P41 support fragment in the local window."""
import csv, pathlib, sys
import pysam

samples_csv, supports_tsv, window, outdir = sys.argv[1:5]
out = pathlib.Path(outdir); out.mkdir(parents=True, exist_ok=True)
chrom, interval = window.split(":"); start, end = map(int, interval.split("-"))
by_sample = {}
for r in csv.DictReader(open(supports_tsv), delimiter="\t"):
    by_sample.setdefault(r["sample"], {}).setdefault(r["qname"], []).append(r["candidate_id"])
manifest = []
for sample in csv.DictReader(open(samples_csv)):
    sid, bam_path = sample["sample_id"], sample["bam_path"]
    wanted = by_sample[sid]
    selected = {}
    with pysam.AlignmentFile(bam_path) as bam:
        for read in bam.fetch(chrom, start - 1, end):
            if read.query_name not in wanted or read.is_secondary or not read.query_sequence: continue
            key = (read.query_name, 1 if read.is_read1 else 2 if read.is_read2 else 0)
            # Prefer a primary record when present in the fetched window; otherwise retain a
            # non-secondary record (including a supplementary record) for this read end.
            if key not in selected or (selected[key].is_supplementary and not read.is_supplementary): selected[key] = read
    with (out / f"{sid}.fastq").open("w") as fq:
        for (qname, mate), read in sorted(selected.items()):
            seq = read.query_sequence
            qual = read.qual or "I" * len(seq)
            for cid in wanted[qname]:
                tag = f"{cid}|{qname}|{mate}"
                fq.write(f"@{tag}\n{seq}\n+\n{qual}\n")
                manifest.append({"candidate_id": cid, "sample": sid, "qname": qname, "mate": mate, "fastq_id": tag})
    print(f"{sid}: extracted {len(selected)} read ends for {len(wanted)} source fragments")
with (out / "rescue_read_manifest.tsv").open("w", newline="") as h:
    w = csv.DictWriter(h, fieldnames=["candidate_id", "sample", "qname", "mate", "fastq_id"], delimiter="\t", lineterminator="\n"); w.writeheader(); w.writerows(manifest)
if not manifest: raise SystemExit("no candidate reads were extracted")
