#!/usr/bin/env python3
"""Recover unmapped mates that cannot be randomly retrieved from a BAM index.

Coordinate-sorted BAM indexes cannot seek to unmapped reads.  For the explicit
QNAMEs where pysam's mate lookup reports `is unmapped`, this script scans only
the affected sample BAMs to EOF, records every matching alignment, and writes
retrieved unmapped sequences as additional single-end local-assembly input.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import json
from collections import defaultdict
from pathlib import Path

import pysam


def read_end(read):
    return "R1" if read.is_read1 else "R2" if read.is_read2 else "U"


def record_class(read):
    if read.is_unmapped:
        return "unmapped"
    if read.is_secondary:
        return "secondary"
    if read.is_supplementary:
        return "supplementary"
    return "primary"


def fastq_write(handle, qname, end, sequence, qualities):
    quality = qualities if qualities and len(qualities) == len(sequence) else "I" * len(sequence)
    handle.write(f"@{qname}|{end}|unmapped_mate_rescue\n{sequence}\n+\n{quality}\n")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--inputs-dir", required=True, type=Path)
    parser.add_argument("--sample-table", required=True, type=Path)
    args = parser.parse_args()
    failures_path = args.inputs_dir / "mate_lookup_failures.tsv"
    failures = list(csv.DictReader(failures_path.open(), delimiter="\t"))
    if not failures:
        raise SystemExit("No mate lookup failures to rescue")
    bam_by_sample = {row["sample_id"]: row["bam_path"]
                     for row in csv.DictReader(args.sample_table.open())}
    requested = defaultdict(set)
    locus_by_sample_qname = {}
    for row in failures:
        requested[row["sample_id"]].add(row["qname"])
        locus_by_sample_qname[(row["sample_id"], row["qname"])] = row["locus_id"]
    rescue_root = args.inputs_dir / "unmapped_mate_rescue"
    rescue_root.mkdir(exist_ok=False)
    ledger, recovered = [], {}
    scan_summary = []
    for sample_id, qnames in sorted(requested.items()):
        bam_path = bam_by_sample[sample_id]
        print(f"START sample={sample_id} target_qnames={len(qnames)} bam={bam_path}", flush=True)
        scanned, matched = 0, 0
        with pysam.AlignmentFile(bam_path, "rb") as bam:
            for read in bam.fetch(until_eof=True):
                scanned += 1
                if scanned % 5_000_000 == 0:
                    print(f"PROGRESS sample={sample_id} scanned_records={scanned} matched_records={matched}", flush=True)
                if read.query_name not in qnames:
                    continue
                matched += 1
                locus_id = locus_by_sample_qname[(sample_id, read.query_name)]
                seq = read.query_sequence or ""
                row = {
                    "locus_id": locus_id, "sample_id": sample_id, "qname": read.query_name,
                    "read_end": read_end(read), "record_class": record_class(read), "flag": read.flag,
                    "chrom": read.reference_name or "*",
                    "start_1based": read.reference_start + 1 if read.reference_start >= 0 else "",
                    "mapq": read.mapping_quality, "cigar": read.cigarstring or "*",
                    "has_sequence": int(bool(seq)), "stored_sequence_length": len(seq),
                }
                ledger.append(row)
                if read.is_unmapped and seq:
                    key = (locus_id, read.query_name, read_end(read))
                    old = recovered.get(key)
                    if old is None or len(seq) > len(old["sequence"]):
                        recovered[key] = {"sequence": seq, "qualities": read.qual or "", **row}
        scan_summary.append({"sample_id": sample_id, "bam_path": bam_path, "target_qnames": len(qnames),
                             "scanned_records": scanned, "matched_records": matched})
        print(f"DONE sample={sample_id} scanned_records={scanned} matched_records={matched}", flush=True)
    ledger_fields = ["locus_id", "sample_id", "qname", "read_end", "record_class", "flag", "chrom",
                     "start_1based", "mapq", "cigar", "has_sequence", "stored_sequence_length"]
    with (rescue_root / "unmapped_mate_recovery_ledger.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=ledger_fields, delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(ledger)
    by_locus = defaultdict(list)
    for payload in recovered.values():
        by_locus[payload["locus_id"]].append(payload)
    for locus_id, records in by_locus.items():
        out = rescue_root / "loci" / locus_id
        out.mkdir(parents=True, exist_ok=True)
        with gzip.open(out / "reads_recovered_mates.fastq.gz", "wt") as handle:
            for payload in sorted(records, key=lambda x: (x["qname"], x["read_end"])):
                fastq_write(handle, payload["qname"], payload["read_end"], payload["sequence"], payload["qualities"])
    status_rows = []
    for failure in failures:
        key_prefix = (failure["locus_id"], failure["qname"])
        hits = [value for key, value in recovered.items() if key[:2] == key_prefix]
        status_rows.append({**failure,
                            "recovery_status": "recovered_unmapped_sequence" if hits else "not_recovered_after_full_bam_scan",
                            "recovered_read_ends": ",".join(sorted({item["read_end"] for item in hits})),
                            "recovered_sequence_lengths": ",".join(str(item["stored_sequence_length"]) for item in hits)})
    status_fields = list(failures[0]) + ["recovery_status", "recovered_read_ends", "recovered_sequence_lengths"]
    with (rescue_root / "unmapped_mate_recovery_status.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=status_fields, delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(status_rows)
    (rescue_root / "audit.json").write_text(json.dumps({
        "requested_unmapped_mates": len(failures), "recovered_unmapped_read_ends": len(recovered),
        "not_recovered_qnames": sum(row["recovery_status"] != "recovered_unmapped_sequence" for row in status_rows),
        "scan_summary": scan_summary,
        "interpretation": "Full BAM scan is used because coordinate indexes cannot seek unmapped reads. A failure to recover after this scan is reported as unrecovered, not SV-negative."
    }, indent=2) + "\n")
    (rescue_root / ".complete").write_text("completed\n")


if __name__ == "__main__":
    main()
