#!/usr/bin/env python3
"""Fragment-aware clustering of ungated P42 junction evidence.

This is candidate generation only, not SV calling.  Every emitted raw
alignment-context row remains traceable to its QNAME, read end, and record
class.  At a cluster, QNAME/template is the primary independent-fragment
support unit; QNAME+read-end is reported separately as read-level support.
Thus primary, secondary, and supplementary records are not discarded or
mistaken for extra fragments, and paired R1/R2 are not counted as two
independent templates.
"""
from __future__ import annotations
import argparse, csv, gzip, json
from collections import defaultdict
from pathlib import Path


def bin_start(pos: int, width: int) -> int:
    return ((pos - 1) // width) * width + 1


def sa_destinations(tag: str):
    for entry in tag.rstrip(";").split(";"):
        fields = entry.split(",")
        if len(fields) >= 4:
            try: yield fields[0], int(fields[1]), fields[2], fields[3]
            except ValueError: pass


def keys(row):
    typ = row["evidence_type"]
    chrom = row["chrom"]
    start = int(row["start_1based"] or 0)
    end = int(row["end_1based"] or 0)
    if typ == "softclip":
        bp = int(row["breakend_1based"])
        yield f"softclip|{chrom}:{bin_start(bp, 100)}|{row['clip_side']}"
    elif typ == "sa_split":
        for other_chrom, other_pos, other_strand, other_cigar in sa_destinations(row["sa_tag"]):
            yield (f"sa_split|{chrom}:{bin_start(start, 100)}:{row['strand']}|"
                   f"{other_chrom}:{bin_start(other_pos, 100)}:{other_strand}")
    elif typ == "supplementary_alignment":
        yield f"supplementary|{chrom}:{bin_start(start, 100)}-{bin_start(end, 100)}:{row['strand']}"
    elif typ == "discordant_pair":
        mate_chrom = row["mate_chrom"] or "*"
        mate_start = int(row["mate_start_1based"] or 0)
        yield f"discordant|{chrom}:{bin_start(start, 1000)}|{mate_chrom}:{bin_start(mate_start, 1000) if mate_start else 0}"


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--input-dir", required=True, type=Path)
    p.add_argument("--out-dir", required=True, type=Path)
    args = p.parse_args(); args.out_dir.mkdir(parents=True, exist_ok=True)
    sample_dirs = sorted(path for path in args.input_dir.iterdir() if path.is_dir())
    if len(sample_dirs) != 18: raise SystemExit(f"Expected 18 evidence directories, found {len(sample_dirs)}")
    raw_rows, clusters = [], defaultdict(
        lambda: {"fragments": set(), "read_ends": set(), "classes": set(), "mapqs": [], "records": 0}
    )
    for sample_dir in sample_dirs:
        evidence = sample_dir / "junction_evidence_records.tsv.gz"
        with gzip.open(evidence, "rt", newline="") as h:
            for row in csv.DictReader(h, delimiter="\t"):
                fragment_id = row["qname"]
                read_end_id = f"{row['qname']}|{row['read_end']}"
                for key in keys(row):
                    raw_rows.append({"sample_id": row["sample_id"], "cluster_key": key,
                                     "fragment_id": fragment_id, "read_end_id": read_end_id,
                                     "evidence_type": row["evidence_type"], "record_class": row["record_class"],
                                     "mapq": row["mapq"], "is_duplicate": row["is_duplicate"],
                                     "qname": row["qname"], "read_end": row["read_end"], "flag": row["flag"]})
                    c = clusters[(row["sample_id"], key)]
                    c["fragments"].add(fragment_id); c["read_ends"].add(read_end_id)
                    c["classes"].add(row["record_class"])
                    c["mapqs"].append(int(row["mapq"])); c["records"] += 1
    raw_fields = ["sample_id", "cluster_key", "fragment_id", "read_end_id", "evidence_type", "record_class",
                  "mapq", "is_duplicate", "qname", "read_end", "flag"]
    with gzip.open(args.out_dir / "support_unit_to_cluster.tsv.gz", "wt", newline="") as h:
        w = csv.DictWriter(h, fieldnames=raw_fields, delimiter="\t", lineterminator="\n"); w.writeheader(); w.writerows(raw_rows)
    summary = []
    for (sample, key), c in sorted(clusters.items()):
        summary.append({"sample_id": sample, "cluster_key": key, "evidence_type": key.split("|", 1)[0],
                        "independent_fragment_count": len(c["fragments"]),
                        "independent_read_end_count": len(c["read_ends"]), "record_count": c["records"],
                        "record_classes": ",".join(sorted(c["classes"])), "max_mapq": max(c["mapqs"]),
                        "median_mapq": sorted(c["mapqs"])[len(c["mapqs"]) // 2]})
    fields = ["sample_id", "cluster_key", "evidence_type", "independent_fragment_count",
              "independent_read_end_count", "record_count", "record_classes", "max_mapq", "median_mapq"]
    with (args.out_dir / "candidate_clusters.tsv").open("w", newline="") as h:
        w = csv.DictWriter(h, fieldnames=fields, delimiter="\t", lineterminator="\n"); w.writeheader(); w.writerows(summary)
    (args.out_dir / "audit.json").write_text(json.dumps({
        "samples": len(sample_dirs),
        "support_unit_cluster_rows": len(raw_rows),
        "candidate_clusters": len(summary),
        "primary_support_unit": "unique (sample_id, QNAME/template): independent fragment",
        "read_level_support_unit": "unique (sample_id, QNAME, read_end): sequenced read end",
        "threshold_rule": "Use independent_fragment_count, not independent_read_end_count, for the primary >=3 independent-fragment threshold.",
        "rule": "100-bp bins for clip/SA anchors, 1-kb bins for discordant pairs; clusters are discovery candidates, not SV calls."
    }, indent=2) + "\n")
    (args.out_dir / ".complete").write_text("completed\n")

if __name__ == "__main__": main()
