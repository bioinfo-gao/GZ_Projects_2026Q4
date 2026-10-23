#!/usr/bin/env python3
"""Prepare auditable per-sample local-assembly inputs for P42 priority bins.

Priority bins are merged only when their same-sample chr11 windows overlap.
Each locus receives all mapped records in that local window plus mate lookups
for the priority-evidence QNAMEs.  Alignment records are never discarded from
the ledger: FASTQ de-duplication is only a sequence-input operation, using the
longest stored sequence per QNAME/read-end and retaining every alternative
record's FLAG, CIGAR, coordinate, MAPQ, and source in the ledger.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import json
from collections import Counter, defaultdict
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


def has_hard_clip(read):
    return int(any(op == 5 for op, _ in (read.cigartuples or [])))


def record_key(read):
    return "|".join(map(str, [read.query_name, read_end(read), read.flag,
                                read.reference_name or "*", read.reference_start,
                                read.reference_end, read.cigarstring or "*"]))


def row_for_record(read, source):
    seq = read.query_sequence or ""
    return {
        "qname": read.query_name,
        "read_end": read_end(read),
        "record_class": record_class(read),
        "flag": read.flag,
        "chrom": read.reference_name or "*",
        "start_1based": read.reference_start + 1 if read.reference_start >= 0 else "",
        "end_1based": read.reference_end if read.reference_end is not None else "",
        "mapq": read.mapping_quality,
        "cigar": read.cigarstring or "*",
        "is_duplicate": int(read.is_duplicate),
        "is_qcfail": int(read.is_qcfail),
        "has_sequence": int(bool(seq)),
        "stored_sequence_length": len(seq),
        "has_hard_clip": has_hard_clip(read),
        "source": source,
    }


def build_loci(manifest, radius):
    by_sample = defaultdict(list)
    for row in manifest:
        if row["anchor_chrom"] != "chr11" or not row["anchor_1based"]:
            raise ValueError(f"Missing/invalid chr11 anchor in manifest row: {row['cluster_key']}")
        pos = int(row["anchor_1based"])
        by_sample[row["sample_id"]].append((pos, row))
    loci, candidate_to_locus = [], {}
    for sample_id, entries in sorted(by_sample.items()):
        current = None
        for pos, row in sorted(entries, key=lambda x: (x[0], x[1]["cluster_key"])):
            start, end = max(1, pos - radius), pos + radius
            if current is None or start > current["end"]:
                if current is not None:
                    loci.append(current)
                current = {"sample_id": sample_id, "chrom": "chr11", "start": start, "end": end,
                           "rows": [row]}
            else:
                current["end"] = max(current["end"], end)
                current["rows"].append(row)
        if current is not None:
            loci.append(current)
    numbered = []
    for index, locus in enumerate(loci, start=1):
        locus_id = f"{locus['sample_id']}_L{index:03d}_{locus['chrom']}_{locus['start']}_{locus['end']}"
        locus["locus_id"] = locus_id
        numbered.append(locus)
        for row in locus["rows"]:
            candidate_to_locus[(row["sample_id"], row["cluster_key"])] = locus_id
    return numbered, candidate_to_locus


def fastq_write(handle, qname, end, sequence, qualities):
    qual = qualities if qualities and len(qualities) == len(sequence) else "I" * len(sequence)
    handle.write(f"@{qname}|{end}\n{sequence}\n+\n{qual}\n")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--support-links", required=True, type=Path)
    parser.add_argument("--sample-table", required=True, type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--radius-bp", type=int, default=500)
    parser.add_argument("--sample-id", help="Safe scoped run for one sample only")
    args = parser.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    with args.manifest.open(newline="") as handle:
        manifest = list(csv.DictReader(handle, delimiter="\t"))
    if args.sample_id:
        manifest = [row for row in manifest if row["sample_id"] == args.sample_id]
    if not manifest:
        raise SystemExit("No manifest rows available for the requested scope")
    loci, candidate_to_locus = build_loci(manifest, args.radius_bp)

    bam_by_sample = {}
    with args.sample_table.open(newline="") as handle:
        for row in csv.DictReader(handle):
            if not args.sample_id or row["sample_id"] == args.sample_id:
                bam_by_sample[row["sample_id"]] = Path(row["bam_path"])
    missing = sorted({locus["sample_id"] for locus in loci} - set(bam_by_sample))
    if missing:
        raise SystemExit(f"No BAM paths for manifest samples: {','.join(missing)}")

    qnames_by_locus = defaultdict(set)
    with gzip.open(args.support_links, "rt", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            locus_id = candidate_to_locus.get((row["sample_id"], row["cluster_key"]))
            if locus_id:
                qnames_by_locus[locus_id].add(row["qname"])

    loci_fields = ["locus_id", "sample_id", "chrom", "start_1based", "end_1based", "width_bp",
                   "priority_cluster_bins", "priority_evidence_qnames"]
    with (args.out_dir / "assembly_loci.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=loci_fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for locus in loci:
            writer.writerow({"locus_id": locus["locus_id"], "sample_id": locus["sample_id"],
                             "chrom": locus["chrom"], "start_1based": locus["start"], "end_1based": locus["end"],
                             "width_bp": locus["end"] - locus["start"] + 1,
                             "priority_cluster_bins": len(locus["rows"]),
                             "priority_evidence_qnames": len(qnames_by_locus[locus["locus_id"]])})

    candidate_fields = list(manifest[0]) + ["locus_id"]
    with (args.out_dir / "candidate_to_locus.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=candidate_fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in manifest:
            row = dict(row)
            row["locus_id"] = candidate_to_locus[(row["sample_id"], row["cluster_key"])]
            writer.writerow(row)

    audit_rows, mate_failures, global_counts = [], [], Counter()
    ledger_fields = ["locus_id", "sample_id", "is_priority_evidence_qname", "qname", "read_end",
                     "record_class", "flag", "chrom", "start_1based", "end_1based", "mapq", "cigar",
                     "is_duplicate", "is_qcfail", "has_sequence", "stored_sequence_length", "has_hard_clip", "source"]
    for locus in loci:
        locus_dir = args.out_dir / "loci" / locus["locus_id"]
        locus_dir.mkdir(parents=True, exist_ok=True)
        wanted_qnames = qnames_by_locus[locus["locus_id"]]
        records = {}
        sequence_by_end = {}
        counts = Counter()
        bam_path = bam_by_sample[locus["sample_id"]]
        with pysam.AlignmentFile(str(bam_path), "rb") as bam:
            primary_for_mate = {}
            for read in bam.fetch(locus["chrom"], locus["start"] - 1, locus["end"]):
                key = record_key(read)
                records[key] = {**row_for_record(read, "context_window"), "sources": {"context_window"}}
                counts["context_records"] += 1
                if read.query_name in wanted_qnames and not read.is_unmapped and not read.is_secondary and not read.is_supplementary:
                    primary_for_mate.setdefault(read.query_name, read)
                if read.query_sequence:
                    end_key = (read.query_name, read_end(read))
                    old = sequence_by_end.get(end_key)
                    if old is None or len(read.query_sequence) > len(old["sequence"]):
                        sequence_by_end[end_key] = {"sequence": read.query_sequence,
                                                    "qualities": read.qual or "",
                                                    "source": "context_window"}
            for qname, read in primary_for_mate.items():
                try:
                    mate = bam.mate(read)
                except (ValueError, OSError) as exc:
                    counts["mate_lookup_failed"] += 1
                    mate_failures.append({"locus_id": locus["locus_id"], "sample_id": locus["sample_id"],
                                          "qname": qname, "source_read_end": read_end(read), "source_flag": read.flag,
                                          "source_chrom": read.reference_name or "*",
                                          "source_start_1based": read.reference_start + 1 if read.reference_start >= 0 else "",
                                          "error_type": type(exc).__name__, "error_message": str(exc)})
                    continue
                key = record_key(mate)
                if key in records:
                    records[key]["sources"].add("mate_lookup_of_priority_qname")
                else:
                    records[key] = {**row_for_record(mate, "mate_lookup_of_priority_qname"),
                                    "sources": {"mate_lookup_of_priority_qname"}}
                counts["mate_lookup_records"] += 1
                if mate.query_sequence:
                    end_key = (mate.query_name, read_end(mate))
                    old = sequence_by_end.get(end_key)
                    if old is None or len(mate.query_sequence) > len(old["sequence"]):
                        sequence_by_end[end_key] = {"sequence": mate.query_sequence,
                                                    "qualities": mate.qual or "",
                                                    "source": "mate_lookup_of_priority_qname"}
        with gzip.open(locus_dir / "alignment_record_ledger.tsv.gz", "wt", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=ledger_fields, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            for record in sorted(records.values(), key=lambda x: (x["qname"], x["read_end"], x["flag"], str(x["start_1based"]))):
                record["locus_id"] = locus["locus_id"]
                record["sample_id"] = locus["sample_id"]
                record["is_priority_evidence_qname"] = int(record["qname"] in wanted_qnames)
                record["source"] = ";".join(sorted(record.pop("sources")))
                writer.writerow(record)
                counts[f"ledger_{record['record_class']}"] += 1
                counts[f"ledger_has_sequence_{record['has_sequence']}"] += 1
                counts[f"ledger_hard_clip_{record['has_hard_clip']}"] += 1
        paired, unpaired = {}, {}
        by_qname = defaultdict(dict)
        for (qname, end), payload in sequence_by_end.items():
            by_qname[qname][end] = payload
        for qname, ends in by_qname.items():
            if "R1" in ends and "R2" in ends:
                paired[qname] = ends
            else:
                unpaired[qname] = ends
        with gzip.open(locus_dir / "reads_R1.fastq.gz", "wt") as r1, \
             gzip.open(locus_dir / "reads_R2.fastq.gz", "wt") as r2, \
             gzip.open(locus_dir / "reads_unpaired.fastq.gz", "wt") as single:
            for qname, ends in sorted(paired.items()):
                fastq_write(r1, qname, "R1", ends["R1"]["sequence"], ends["R1"]["qualities"])
                fastq_write(r2, qname, "R2", ends["R2"]["sequence"], ends["R2"]["qualities"])
            for qname, ends in sorted(unpaired.items()):
                for end, payload in sorted(ends.items()):
                    fastq_write(single, qname, end, payload["sequence"], payload["qualities"])
        counts["fastq_paired_templates"] = len(paired)
        counts["fastq_unpaired_read_ends"] = sum(len(x) for x in unpaired.values())
        counts["priority_evidence_qnames"] = len(wanted_qnames)
        audit_rows.append({"locus_id": locus["locus_id"], "sample_id": locus["sample_id"],
                           "chrom": locus["chrom"], "start_1based": locus["start"], "end_1based": locus["end"],
                           **dict(sorted(counts.items()))})
        global_counts.update(counts)
    fixed = ["locus_id", "sample_id", "chrom", "start_1based", "end_1based"]
    all_fields = fixed + sorted({key for row in audit_rows for key in row} - set(fixed))
    with (args.out_dir / "assembly_input_audit.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=all_fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(audit_rows)
    with (args.out_dir / "mate_lookup_failures.tsv").open("w", newline="") as handle:
        fields = ["locus_id", "sample_id", "qname", "source_read_end", "source_flag", "source_chrom",
                  "source_start_1based", "error_type", "error_message"]
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(mate_failures)
    (args.out_dir / "audit.json").write_text(json.dumps({
        "scope_sample_id": args.sample_id or "all_18_samples",
        "radius_bp": args.radius_bp,
        "n_loci": len(loci),
        "mate_lookup_failure_rows": len(mate_failures),
        "primary_support_unit": "unique QNAME/template; no alignment record removed",
        "fastq_sequence_rule": "Within a locus, one FASTQ sequence per QNAME/read-end: select the longest nonempty stored BAM SEQ. Ledger retains all original records; hard clipping is audited and nonempty SEQ is not assumed complete.",
        "global_counts": dict(sorted(global_counts.items())),
    }, indent=2) + "\n")
    (args.out_dir / ".complete").write_text("completed\n")


if __name__ == "__main__":
    main()
