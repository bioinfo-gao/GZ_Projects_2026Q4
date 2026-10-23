#!/usr/bin/env python3
"""Extract all P42 junction-relevant BAM evidence without candidate/QNAME gating."""
from __future__ import annotations

import argparse
import csv
import gzip
import json
from collections import Counter
from pathlib import Path

import pysam

FIELDS = [
    "sample_id", "evidence_type", "qname", "read_end", "record_class", "flag",
    "is_duplicate", "is_qcfail", "mapq", "chrom", "start_1based", "end_1based",
    "strand", "cigar", "has_sequence", "query_length", "clip_side", "clip_length",
    "breakend_1based", "mate_chrom", "mate_start_1based", "template_length", "sa_tag",
]


def read_end(read):
    return "R1" if read.is_read1 else "R2" if read.is_read2 else "U"


def alignment_class(read):
    if read.is_unmapped: return "unmapped"
    if read.is_secondary: return "secondary"
    if read.is_supplementary: return "supplementary"
    return "primary"


def clips(read):
    c = read.cigartuples or []
    return (c[0][1] if c and c[0][0] == 4 else 0,
            c[-1][1] if c and c[-1][0] == 4 else 0)


def base_row(sample_id, read):
    left, right = clips(read)
    return {
        "sample_id": sample_id, "qname": read.query_name, "read_end": read_end(read),
        "record_class": alignment_class(read), "flag": read.flag,
        "is_duplicate": int(read.is_duplicate), "is_qcfail": int(read.is_qcfail),
        "mapq": read.mapping_quality, "chrom": read.reference_name or "*",
        "start_1based": read.reference_start + 1 if read.reference_start >= 0 else "",
        "end_1based": read.reference_end if read.reference_end is not None else "",
        "strand": "-" if read.is_reverse else "+", "cigar": read.cigarstring or "*",
        "has_sequence": int(read.query_sequence is not None), "query_length": read.query_length or 0,
        "clip_side": "", "clip_length": "", "breakend_1based": "",
        "mate_chrom": read.next_reference_name or "*",
        "mate_start_1based": read.next_reference_start + 1 if read.next_reference_start >= 0 else "",
        "template_length": read.template_length,
        "sa_tag": read.get_tag("SA") if read.has_tag("SA") else "",
        "_left_clip": left, "_right_clip": right,
    }


def evidence_rows(sample_id, read):
    row = base_row(sample_id, read)
    rows = []
    # 20 nt is an evidence-extraction threshold, not a proof threshold.  All
    # records with shorter clips remain traceable in the stage-01 record ledger.
    for side, length, endpoint in (("left", row["_left_clip"], row["start_1based"]),
                                   ("right", row["_right_clip"], row["end_1based"])):
        if length >= 20:
            clip = {key: value for key, value in row.items() if not key.startswith("_")}
            clip.update(evidence_type="softclip", clip_side=side, clip_length=length,
                        breakend_1based=endpoint)
            rows.append(clip)
    if row["sa_tag"]:
        sa = {key: value for key, value in row.items() if not key.startswith("_")}
        sa.update(evidence_type="sa_split")
        rows.append(sa)
    if read.is_supplementary:
        supp = {key: value for key, value in row.items() if not key.startswith("_")}
        supp.update(evidence_type="supplementary_alignment")
        rows.append(supp)
    discordant = (read.mate_is_unmapped or
                  (read.next_reference_name and read.next_reference_name != read.reference_name) or
                  (not read.mate_is_unmapped and read.next_reference_start >= 0 and
                   abs(read.next_reference_start - read.reference_start) > 1000) or
                  not read.is_proper_pair)
    if discordant:
        pair = {key: value for key, value in row.items() if not key.startswith("_")}
        pair.update(evidence_type="discordant_pair")
        rows.append(pair)
    return rows


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--sample-id", required=True)
    p.add_argument("--bam", required=True, type=Path)
    p.add_argument("--out-dir", required=True, type=Path)
    p.add_argument("--chrom", default="chr11")
    p.add_argument("--start", type=int, default=114554000)
    p.add_argument("--end", type=int, default=114675500)
    args = p.parse_args(); args.out_dir.mkdir(parents=True, exist_ok=True)
    counts = Counter(); output = args.out_dir / "junction_evidence_records.tsv.gz"
    with pysam.AlignmentFile(str(args.bam), "rb") as bam, gzip.open(output, "wt", newline="") as h:
        writer = csv.DictWriter(h, fieldnames=FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for read in bam.fetch(args.chrom, args.start - 1, args.end):
            counts["records_in_universe"] += 1
            counts[f"records_{alignment_class(read)}"] += 1
            for row in evidence_rows(args.sample_id, read):
                writer.writerow(row); counts["evidence_rows_total"] += 1
                counts[f"evidence_{row['evidence_type']}"] += 1
                counts[f"evidence_{row['record_class']}"] += 1
    audit = {"sample_id": args.sample_id, "region": f"{args.chrom}:{args.start}-{args.end}",
             "filter_note": "No record class, MAPQ, duplicate or QNAME gate was applied. Clip >=20 nt only controls which existing records are emitted as softclip evidence; all records remain in stage-01 ledger.",
             "counts": dict(sorted(counts.items()))}
    (args.out_dir / "audit.json").write_text(json.dumps(audit, indent=2, sort_keys=True) + "\n")
    (args.out_dir / ".complete").write_text("completed\n")


if __name__ == "__main__":
    main()
