#!/usr/bin/env python3
"""Write all and summary placement ledgers for P42 assembled contigs.

No placement is declared true merely because it is primary.  The output is a
mapping-evidence ledger used by later breakpoint/ambiguity classification.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import json
from collections import Counter, defaultdict
from pathlib import Path

import pysam


def fasta_lengths(path):
    lengths, name, chunks = {}, None, []
    with path.open() as handle:
        for line in handle:
            line = line.rstrip()
            if line.startswith(">"):
                if name is not None:
                    lengths[name] = sum(map(len, chunks))
                name, chunks = line[1:].split()[0], []
            else:
                chunks.append(line)
    if name is not None:
        lengths[name] = sum(map(len, chunks))
    return lengths


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--contigs", required=True, type=Path)
    parser.add_argument("--sam", required=True, type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    args = parser.parse_args()
    contig_lengths = fasta_lengths(args.contigs)
    ledger_fields = ["contig_id", "contig_length", "alignment_class", "flag", "chrom", "start_1based",
                     "end_1based", "strand", "mapq", "cigar", "alignment_score", "edit_distance"]
    rows, grouped, counts = [], defaultdict(list), Counter()
    with pysam.AlignmentFile(str(args.sam), "r") as alignment:
        for read in alignment:
            if read.is_unmapped:
                alignment_class, chrom, start, end = "unmapped", "*", "", ""
            elif read.is_secondary:
                alignment_class, chrom = "secondary", read.reference_name
                start, end = read.reference_start + 1, read.reference_end
            elif read.is_supplementary:
                alignment_class, chrom = "supplementary", read.reference_name
                start, end = read.reference_start + 1, read.reference_end
            else:
                alignment_class, chrom = "primary", read.reference_name
                start, end = read.reference_start + 1, read.reference_end
            row = {
                "contig_id": read.query_name, "contig_length": contig_lengths.get(read.query_name, ""),
                "alignment_class": alignment_class, "flag": read.flag, "chrom": chrom,
                "start_1based": start, "end_1based": end,
                "strand": "-" if read.is_reverse else "+", "mapq": read.mapping_quality,
                "cigar": read.cigarstring or "*",
                "alignment_score": read.get_tag("AS") if read.has_tag("AS") else "",
                "edit_distance": read.get_tag("NM") if read.has_tag("NM") else "",
            }
            rows.append(row); grouped[read.query_name].append(row)
            counts[f"alignment_{alignment_class}"] += 1
    with gzip.open(args.out_dir / "contig_all_placement_ledger.tsv.gz", "wt", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=ledger_fields, delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)
    summary_fields = ["contig_id", "contig_length", "all_alignment_records", "mapped_alignment_records",
                      "primary_records", "secondary_records", "supplementary_records", "distinct_loci",
                      "primary_chrom", "primary_start_1based", "primary_end_1based", "primary_strand",
                      "primary_mapq", "primary_alignment_score", "alternative_loci"]
    summaries = []
    for contig_id, group in sorted(grouped.items()):
        mapped = [row for row in group if row["alignment_class"] != "unmapped"]
        primary = next((row for row in group if row["alignment_class"] == "primary"), {})
        loci = sorted({f"{row['chrom']}:{row['start_1based']}-{row['end_1based']}:{row['strand']}"
                       for row in mapped})
        alternatives = [loc for loc in loci if loc != f"{primary.get('chrom','')}:{primary.get('start_1based','')}-{primary.get('end_1based','')}:{primary.get('strand','')}"]
        summaries.append({
            "contig_id": contig_id, "contig_length": contig_lengths.get(contig_id, ""),
            "all_alignment_records": len(group), "mapped_alignment_records": len(mapped),
            "primary_records": sum(r["alignment_class"] == "primary" for r in group),
            "secondary_records": sum(r["alignment_class"] == "secondary" for r in group),
            "supplementary_records": sum(r["alignment_class"] == "supplementary" for r in group),
            "distinct_loci": len(loci), "primary_chrom": primary.get("chrom", ""),
            "primary_start_1based": primary.get("start_1based", ""), "primary_end_1based": primary.get("end_1based", ""),
            "primary_strand": primary.get("strand", ""), "primary_mapq": primary.get("mapq", ""),
            "primary_alignment_score": primary.get("alignment_score", ""),
            "alternative_loci": ";".join(alternatives),
        })
    with (args.out_dir / "contig_placement_summary.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=summary_fields, delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(summaries)
    (args.out_dir / "audit.json").write_text(json.dumps({
        "input_contigs": len(contig_lengths), "contigs_with_sam_records": len(grouped),
        "alignment_record_counts": dict(sorted(counts.items())),
        "interpretation_boundary": "Primary/secondary/supplementary labels are retained mapping hypotheses. This file does not designate a true locus or SV breakpoint."
    }, indent=2) + "\n")


if __name__ == "__main__":
    main()
