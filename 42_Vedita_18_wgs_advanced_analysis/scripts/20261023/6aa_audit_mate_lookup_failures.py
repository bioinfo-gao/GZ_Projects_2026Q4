#!/usr/bin/env python3
"""Create a QNAME-level mate-lookup failure audit for an existing P42 6a run."""
from __future__ import annotations

import argparse
import csv
import gzip
import json
from collections import defaultdict
from pathlib import Path

import pysam


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--inputs-dir", required=True, type=Path)
    parser.add_argument("--support-links", required=True, type=Path)
    parser.add_argument("--sample-table", required=True, type=Path)
    args = parser.parse_args()
    loci = list(csv.DictReader((args.inputs_dir / "assembly_loci.tsv").open(), delimiter="\t"))
    candidate_rows = csv.DictReader((args.inputs_dir / "candidate_to_locus.tsv").open(), delimiter="\t")
    cluster_to_locus = {(row["sample_id"], row["cluster_key"]): row["locus_id"] for row in candidate_rows}
    wanted_by_locus = defaultdict(set)
    with gzip.open(args.support_links, "rt", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            locus_id = cluster_to_locus.get((row["sample_id"], row["cluster_key"]))
            if locus_id:
                wanted_by_locus[locus_id].add(row["qname"])
    bam_by_sample = {row["sample_id"]: row["bam_path"]
                     for row in csv.DictReader(args.sample_table.open())}
    failures, attempted = [], 0
    for locus in loci:
        sample_id, locus_id = locus["sample_id"], locus["locus_id"]
        wanted = wanted_by_locus[locus_id]
        with pysam.AlignmentFile(bam_by_sample[sample_id], "rb") as bam:
            primary = {}
            for read in bam.fetch(locus["chrom"], int(locus["start_1based"]) - 1, int(locus["end_1based"])):
                if read.query_name in wanted and not read.is_unmapped and not read.is_secondary and not read.is_supplementary:
                    primary.setdefault(read.query_name, read)
            for qname, read in primary.items():
                attempted += 1
                try:
                    bam.mate(read)
                except (ValueError, OSError) as exc:
                    failures.append({
                        "locus_id": locus_id, "sample_id": sample_id, "qname": qname,
                        "source_flag": read.flag, "source_chrom": read.reference_name or "*",
                        "source_start_1based": read.reference_start + 1 if read.reference_start >= 0 else "",
                        "error_type": type(exc).__name__, "error_message": str(exc),
                    })
    fields = ["locus_id", "sample_id", "qname", "source_flag", "source_chrom", "source_start_1based",
              "error_type", "error_message"]
    with (args.inputs_dir / "mate_lookup_failures.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(failures)
    (args.inputs_dir / "mate_lookup_failure_audit.json").write_text(json.dumps({
        "attempted_primary_qname_mate_lookups": attempted,
        "failure_rows": len(failures),
        "interpretation": "A mate lookup failure means no retrievable mate record through the BAM index for this attempt. It is not an SV-negative finding and does not erase the source record or its sequence."
    }, indent=2) + "\n")


if __name__ == "__main__":
    main()
