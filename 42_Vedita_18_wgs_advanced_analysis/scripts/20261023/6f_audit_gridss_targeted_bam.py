#!/usr/bin/env python3
"""Create and audit caller-compatible GRIDSS BAMs without erasing raw context.

The GRIDSS ROI extractor intentionally retains every alignment for a selected
fragment.  This script preserves that raw BAM, then excludes only secondary
alignments (FLAG 0x100) from a separate coordinate-sorted caller input.  The
TSV records the impact at record, QNAME/read-end, and template levels.
"""
from __future__ import annotations

import argparse
import csv
from collections import Counter, defaultdict
from pathlib import Path

import pysam


def read_end(read: pysam.AlignedSegment) -> str:
    return "R1" if read.is_read1 else "R2" if read.is_read2 else "U"


def record_class(read: pysam.AlignedSegment) -> str:
    if read.is_unmapped:
        return "unmapped"
    if read.is_secondary:
        return "secondary"
    if read.is_supplementary:
        return "supplementary"
    return "primary"


def metrics(records: list[pysam.AlignedSegment]) -> dict[str, int]:
    qnames = {r.query_name for r in records}
    read_ends = {(r.query_name, read_end(r)) for r in records}
    return {
        "records": len(records),
        "qnames_templates": len(qnames),
        "qname_read_ends": len(read_ends),
        "sequence_available_records": sum(bool(r.query_sequence) for r in records),
        "primary_records": sum(record_class(r) == "primary" for r in records),
        "secondary_records": sum(record_class(r) == "secondary" for r in records),
        "supplementary_records": sum(record_class(r) == "supplementary" for r in records),
        "unmapped_records": sum(record_class(r) == "unmapped" for r in records),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--raw-bam", type=Path, required=True)
    parser.add_argument("--caller-bam", type=Path, required=True)
    parser.add_argument("--audit-tsv", type=Path, required=True)
    parser.add_argument("--threads", type=int, default=2)
    args = parser.parse_args()

    if not args.raw_bam.is_file() or args.raw_bam.stat().st_size == 0:
        raise SystemExit(f"Missing raw extractor BAM: {args.raw_bam}")
    args.caller_bam.parent.mkdir(parents=True, exist_ok=True)
    args.audit_tsv.parent.mkdir(parents=True, exist_ok=True)
    unsorted = args.caller_bam.with_suffix(".unsorted.bam")
    if unsorted.exists() or args.caller_bam.exists():
        raise SystemExit("Caller-compatible output exists; refusing to overwrite historical output")

    raw_records: list[pysam.AlignedSegment] = []
    kept_records: list[pysam.AlignedSegment] = []
    removed_records: list[pysam.AlignedSegment] = []
    with pysam.AlignmentFile(str(args.raw_bam), "rb") as source:
        with pysam.AlignmentFile(str(unsorted), "wb", template=source) as dest:
            for record in source.fetch(until_eof=True):
                raw_records.append(record)
                if record.is_secondary:
                    removed_records.append(record)
                else:
                    kept_records.append(record)
                    dest.write(record)

    pysam.sort("-@", str(args.threads), "-o", str(args.caller_bam), str(unsorted))
    unsorted.unlink()
    pysam.index("-@", str(args.threads), str(args.caller_bam))
    pysam.quickcheck(str(args.caller_bam))

    rows = []
    for category, records, reason in (
        ("raw_gridss_extractor_output", raw_records, "all records for ROI-overlap fragment QNAMEs retained"),
        ("excluded_from_gridss_input", removed_records, "secondary alignment FLAG 0x100 excluded only from caller-compatible BAM"),
        ("gridss_caller_compatible_input", kept_records, "primary, supplementary and unmapped records retained"),
    ):
        rows.append({"sample_id": args.sample_id, "category": category, "reason": reason, **metrics(records)})

    fields = ["sample_id", "category", "reason", "records", "qnames_templates", "qname_read_ends",
              "sequence_available_records", "primary_records", "secondary_records", "supplementary_records",
              "unmapped_records"]
    with args.audit_tsv.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
