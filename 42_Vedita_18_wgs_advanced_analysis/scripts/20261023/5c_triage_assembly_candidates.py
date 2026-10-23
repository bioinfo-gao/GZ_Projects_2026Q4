#!/usr/bin/env python3
"""Create an auditable local-assembly manifest from P42 evidence bins.

This script does not call SVs and does not delete or suppress low-support
evidence.  It marks a bounded set of high-information records for the first
per-sample local-assembly pass, while writing every input cluster to an audit
table with an explicit reason for its tier.  The primary support threshold is
unique QNAME/template (independent_fragment_count); read-end support remains a
separate descriptive field.
"""
from __future__ import annotations

import argparse
import csv
import json
from collections import Counter
from pathlib import Path


def chr_pos(text: str):
    """Return chrom, 1-based bin coordinate from `chrom:position[:strand]`."""
    fields = text.split(":")
    return fields[0], int(fields[1])


def anchor_and_partner(cluster_key: str):
    """Return chr11-side anchor and optional partner; None if unparseable."""
    fields = cluster_key.split("|")
    evidence_type = fields[0]
    if evidence_type in {"sa_split", "discordant"} and len(fields) == 3:
        chrom, pos = chr_pos(fields[1])
        partner_chrom, partner_pos = chr_pos(fields[2])
        return chrom, pos, partner_chrom, partner_pos
    if evidence_type == "softclip" and len(fields) == 3:
        chrom, pos = chr_pos(fields[1])
        return chrom, pos, "", ""
    if evidence_type == "supplementary" and len(fields) == 2:
        chrom, interval_and_strand = fields[1].split(":", 1)
        interval = interval_and_strand.rsplit(":", 1)[0]
        return chrom, int(interval.split("-", 1)[0]), "", ""
    return "", "", "", ""


def caller_hits(sample_id, chrom, pos, callers, tolerance):
    hits = []
    for caller in callers.get(sample_id, []):
        if caller["chrom"] == chrom and abs(caller["pos"] - pos) <= tolerance:
            hits.append(f"{caller['caller']}:{caller['record_id']}")
        elif caller["end_chrom"] == chrom and abs(caller["end_pos"] - pos) <= tolerance:
            hits.append(f"{caller['caller']}:{caller['record_id']}")
    return sorted(set(hits))


def is_distant_or_interchromosomal(chrom, pos, partner_chrom, partner_pos):
    return partner_chrom != chrom or abs(partner_pos - pos) > 5_000


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--clusters", required=True, type=Path)
    parser.add_argument("--caller-records", required=True, type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    parser.add_argument("--caller-tolerance-bp", type=int, default=1_000)
    args = parser.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)

    callers = {}
    with args.caller_records.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["source_readable"] != "1" or row["in_discovery_window"] != "1":
                continue
            callers.setdefault(row["sample_id"], []).append({
                "caller": row["caller"], "record_id": row["record_id"],
                "chrom": row["chrom"], "pos": int(row["pos_1based"]),
                "end_chrom": row["end_chrom"], "end_pos": int(row["end_1based"]),
            })

    audit_rows, assembly_rows = [], []
    counts = Counter()
    with args.clusters.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            n_fragment = int(row["independent_fragment_count"])
            chrom, pos, partner_chrom, partner_pos = anchor_and_partner(row["cluster_key"])
            hit_list = caller_hits(row["sample_id"], chrom, pos, callers, args.caller_tolerance_bp) if chrom else []
            evidence_type = row["evidence_type"]
            explicit_high_support = (
                evidence_type == "sa_split" and n_fragment >= 3
            )
            distant_discordant_high_support = (
                evidence_type == "discordant" and n_fragment >= 3
                and is_distant_or_interchromosomal(chrom, pos, partner_chrom, partner_pos)
            )
            supplementary_high_support = evidence_type == "supplementary" and n_fragment >= 3
            caller_anchored = bool(hit_list) and n_fragment >= 2
            reasons = []
            if explicit_high_support:
                reasons.append("SA_split_with_at_least_3_fragments")
            if distant_discordant_high_support:
                reasons.append("distant_or_interchromosomal_discordant_with_at_least_3_fragments")
            if supplementary_high_support:
                reasons.append("supplementary_with_at_least_3_fragments")
            if caller_anchored:
                reasons.append("within_1kb_of_readable_caller_endpoint_and_at_least_2_fragments")
            if reasons:
                tier = "assembly_priority"
            elif evidence_type in {"sa_split", "discordant", "supplementary"}:
                tier = "explicit_junction_context_not_yet_assembled"
            else:
                tier = "softclip_context_not_yet_assembled"
            enriched = dict(row)
            enriched.update({
                "anchor_chrom": chrom, "anchor_1based": pos,
                "partner_chrom": partner_chrom, "partner_1based": partner_pos,
                "caller_hits_within_1kb": ",".join(hit_list),
                "selection_tier": tier,
                "selection_reason": ";".join(reasons) if reasons else "retained_in_audit_not_selected_for_first_assembly_pass",
            })
            audit_rows.append(enriched)
            if tier == "assembly_priority":
                assembly_rows.append(enriched)
            counts["input_clusters"] += 1
            counts[f"input_{evidence_type}"] += 1
            counts[f"tier_{tier}"] += 1
            counts[f"tier_{tier}_{evidence_type}"] += 1

    fields = [
        "sample_id", "cluster_key", "evidence_type", "independent_fragment_count",
        "independent_read_end_count", "record_count", "record_classes", "max_mapq", "median_mapq",
        "anchor_chrom", "anchor_1based", "partner_chrom", "partner_1based",
        "caller_hits_within_1kb", "selection_tier", "selection_reason",
    ]
    for filename, rows in (("cluster_candidate_audit.tsv", audit_rows),
                           ("assembly_candidate_manifest.tsv", assembly_rows)):
        with (args.out_dir / filename).open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
            writer.writeheader()
            writer.writerows(rows)
    (args.out_dir / "audit.json").write_text(json.dumps({
        "primary_support_unit": "unique (sample_id, QNAME/template): independent fragment",
        "read_level_support_unit": "unique (sample_id, QNAME, read_end): sequenced read end",
        "selection_policy": {
            "assembly_priority": [
                "SA split cluster with >=3 independent fragments",
                "distant/interchromosomal discordant cluster with >=3 independent fragments",
                "supplementary cluster with >=3 independent fragments",
                "any evidence cluster within 1 kb of a readable caller endpoint with >=2 independent fragments",
            ],
            "not_selected": "All remaining evidence clusters remain in cluster_candidate_audit.tsv; not selected for this first bounded assembly pass is not a negative SV finding."
        },
        "counts": dict(sorted(counts.items())),
    }, indent=2) + "\n")
    (args.out_dir / ".complete").write_text("completed\n")


if __name__ == "__main__":
    main()
