#!/usr/bin/env python3
"""Inventory every record in the P42 discovery window without candidate gating.

Writes record-level and molecule-level TSV.GZ ledgers.  It intentionally retains
primary, secondary, supplementary, MAPQ0 and duplicate records.  A molecule is
identified by QNAME plus read-end; record count is never used as molecule support.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import json
import sys
from collections import Counter, defaultdict
from pathlib import Path

import pysam


RECORD_FIELDS = [
    "sample_id", "source", "qname", "read_end", "record_class", "flag",
    "chrom", "start_1based", "end_1based", "strand", "mapq", "cigar",
    "is_duplicate", "is_qcfail", "is_unmapped", "mate_unmapped",
    "mate_chrom", "mate_start_1based", "has_sequence", "query_length",
    "softclip_left", "softclip_right", "has_sa_tag",
]
MOLECULE_FIELDS = [
    "sample_id", "molecule_id", "qname", "read_end", "record_count",
    "primary_count", "secondary_count", "supplementary_count", "mapq0_count",
    "duplicate_record_count", "qcfail_record_count", "sequence_record_count",
    "softclip_record_count", "sa_record_count", "unmapped_mate_rescued_count",
    "max_mapq", "mapped_chromosomes", "mapped_start_min_1based",
    "mapped_end_max_1based",
]


def read_end(read: pysam.AlignedSegment) -> str:
    if read.is_read1:
        return "R1"
    if read.is_read2:
        return "R2"
    return "U"


def record_class(read: pysam.AlignedSegment) -> str:
    if read.is_unmapped:
        return "unmapped"
    if read.is_secondary:
        return "secondary"
    if read.is_supplementary:
        return "supplementary"
    return "primary"


def softclips(cigar: list[tuple[int, int]] | None) -> tuple[int, int]:
    if not cigar:
        return 0, 0
    left = cigar[0][1] if cigar[0][0] == 4 else 0
    right = cigar[-1][1] if cigar[-1][0] == 4 else 0
    return left, right


def row_for_read(sample_id: str, read: pysam.AlignedSegment, source: str) -> dict[str, object]:
    left_clip, right_clip = softclips(read.cigartuples)
    chrom = "*" if read.reference_id < 0 else read.reference_name
    mate_chrom = "*" if read.next_reference_id < 0 else read.next_reference_name
    return {
        "sample_id": sample_id, "source": source, "qname": read.query_name,
        "read_end": read_end(read), "record_class": record_class(read),
        "flag": read.flag, "chrom": chrom,
        "start_1based": "" if read.reference_start < 0 else read.reference_start + 1,
        "end_1based": "" if read.reference_end is None else read.reference_end,
        "strand": "-" if read.is_reverse else "+", "mapq": read.mapping_quality,
        "cigar": read.cigarstring or "*", "is_duplicate": int(read.is_duplicate),
        "is_qcfail": int(read.is_qcfail), "is_unmapped": int(read.is_unmapped),
        "mate_unmapped": int(read.mate_is_unmapped), "mate_chrom": mate_chrom,
        "mate_start_1based": "" if read.next_reference_start < 0 else read.next_reference_start + 1,
        "has_sequence": int(read.query_sequence is not None),
        "query_length": read.query_length or 0, "softclip_left": left_clip,
        "softclip_right": right_clip, "has_sa_tag": int(read.has_tag("SA")),
    }


def update_molecule(molecules: dict, row: dict[str, object]) -> None:
    key = (row["qname"], row["read_end"])
    m = molecules[key]
    m["record_count"] += 1
    m[f"{row['record_class']}_count"] += 1
    m["mapq0_count"] += int(row["mapq"] == 0)
    m["duplicate_record_count"] += int(row["is_duplicate"])
    m["qcfail_record_count"] += int(row["is_qcfail"])
    m["sequence_record_count"] += int(row["has_sequence"])
    m["softclip_record_count"] += int(row["softclip_left"] > 0 or row["softclip_right"] > 0)
    m["sa_record_count"] += int(row["has_sa_tag"])
    m["unmapped_mate_rescued_count"] += int(row["source"] == "mate_rescue")
    m["max_mapq"] = max(m["max_mapq"], int(row["mapq"]))
    if row["chrom"] != "*":
        m["chromosomes"].add(str(row["chrom"]))
    if row["start_1based"] != "":
        m["starts"].append(int(row["start_1based"]))
    if row["end_1based"] != "":
        m["ends"].append(int(row["end_1based"]))


def empty_molecule() -> dict:
    return {
        "record_count": 0, "primary_count": 0, "secondary_count": 0,
        "supplementary_count": 0, "unmapped_count": 0, "mapq0_count": 0,
        "duplicate_record_count": 0, "qcfail_record_count": 0,
        "sequence_record_count": 0, "softclip_record_count": 0,
        "sa_record_count": 0, "unmapped_mate_rescued_count": 0, "max_mapq": 0,
        "chromosomes": set(), "starts": [], "ends": [],
    }


def write_tsv_gz(path: Path, fields: list[str], rows) -> None:
    with gzip.open(path, "wt", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def inventory(sample_id: str, bam_path: Path, chrom: str, start: int, end: int, output: Path) -> None:
    output.mkdir(parents=True, exist_ok=True)
    records_path = output / "records.tsv.gz"
    molecules: dict = defaultdict(empty_molecule)
    counts = Counter()
    rescue_failures = 0
    attempted_rescues: set[tuple[str, str]] = set()

    with pysam.AlignmentFile(str(bam_path), "rb") as bam, gzip.open(records_path, "wt", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=RECORD_FIELDS, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for read in bam.fetch(chrom, start - 1, end):
            row = row_for_read(sample_id, read, "window_fetch")
            writer.writerow(row)
            update_molecule(molecules, row)
            counts["records_total"] += 1
            counts[f"records_{row['record_class']}"] += 1
            counts["records_mapq0"] += int(read.mapping_quality == 0)
            counts["records_has_sequence"] += int(read.query_sequence is not None)
            counts["records_softclipped"] += int(row["softclip_left"] > 0 or row["softclip_right"] > 0)
            counts["records_sa"] += int(read.has_tag("SA"))
            counts["records_duplicate"] += int(read.is_duplicate)
            counts["records_qcfail"] += int(read.is_qcfail)
            if read.mate_is_unmapped and not read.is_unmapped:
                mkey = (read.query_name, read_end(read))
                if mkey not in attempted_rescues:
                    attempted_rescues.add(mkey)
                    counts["unmapped_mate_rescue_attempts"] += 1
                    try:
                        mate = bam.mate(read)
                    except (ValueError, OSError):
                        rescue_failures += 1
                    else:
                        mate_row = row_for_read(sample_id, mate, "mate_rescue")
                        writer.writerow(mate_row)
                        update_molecule(molecules, mate_row)
                        counts["unmapped_mate_rescue_retrieved"] += 1

    molecule_rows = []
    for (qname, end_name), m in sorted(molecules.items()):
        molecule_rows.append({
            "sample_id": sample_id, "molecule_id": f"{qname}|{end_name}",
            "qname": qname, "read_end": end_name,
            **{key: m[key] for key in MOLECULE_FIELDS if key in m},
            "mapped_chromosomes": ",".join(sorted(m["chromosomes"])),
            "mapped_start_min_1based": min(m["starts"]) if m["starts"] else "",
            "mapped_end_max_1based": max(m["ends"]) if m["ends"] else "",
        })
    write_tsv_gz(output / "molecules.tsv.gz", MOLECULE_FIELDS, molecule_rows)
    counts["molecules_total"] = len(molecules)
    counts["unmapped_mate_rescue_failures"] = rescue_failures
    counts["records_with_seq_star"] = counts["records_total"] - counts["records_has_sequence"]
    summary = {
        "sample_id": sample_id, "bam": str(bam_path), "region": f"{chrom}:{start}-{end}",
        "record_counts": dict(sorted(counts.items())),
        "molecule_count": len(molecules),
        "unmapped_mate_rescue_note": "mate() retrieval failures are reported, never treated as absent mate evidence.",
    }
    (output / "audit.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    (output / ".complete").write_text("completed\n")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample-table", required=True, type=Path)
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--chrom", default="chr11")
    parser.add_argument("--start", default=114554000, type=int)
    parser.add_argument("--end", default=114675500, type=int)
    parser.add_argument("--out-root", required=True, type=Path)
    args = parser.parse_args()
    with args.sample_table.open(newline="") as handle:
        samples = {row["sample_id"]: row for row in csv.DictReader(handle)}
    if args.sample_id not in samples:
        raise SystemExit(f"sample absent from canonical table: {args.sample_id}")
    bam = Path(samples[args.sample_id]["bam_path"])
    if not bam.is_file() or not Path(str(bam) + ".bai").is_file():
        raise SystemExit(f"BAM or BAI unavailable: {bam}")
    inventory(args.sample_id, bam, args.chrom, args.start, args.end, args.out_root / args.sample_id)


if __name__ == "__main__":
    main()
