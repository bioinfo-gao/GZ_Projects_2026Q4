#!/usr/bin/env python3
"""Audit and extract P41 Manta/DELLY calls without assuming missing = negative."""
from __future__ import annotations

import argparse
import csv
import re
from pathlib import Path

import pysam

MATE_RE = re.compile(r"[\[\]]([^:\[\]]+):(\d+)[\[\]]")
FIELDS = ["sample_id", "caller", "source_stage", "source_path", "source_readable",
          "filter", "chrom", "pos_1based", "end_chrom", "end_1based", "svtype",
          "svlen", "record_id", "qual", "alt", "in_discovery_window"]


def endpoint(record):
    alt = str(record.alts[0]) if record.alts else "."
    # Manta does not declare CHR2 in its header.  Calling info.get("CHR2")
    # on an undeclared key raises ValueError in pysam rather than returning None.
    if "CHR2" in record.info and record.info["CHR2"]:
        return str(record.info["CHR2"]), int(record.info.get("END", record.pos))
    match = MATE_RE.search(alt)
    if match:
        return match.group(1), int(match.group(2))
    return record.chrom, int(record.stop or record.pos)


def in_window(chrom, pos, end_chrom, end, target_chrom, start, stop):
    return ((chrom == target_chrom and start <= pos <= stop) or
            (end_chrom == target_chrom and start <= end <= stop))


def candidate_files(p41: Path, sid: str):
    return [
        ("Manta", "filtered_result", p41 / "analysis/sv/manta_gw" / sid / "results/variants/diploidSV.vcf.gz"),
        ("Manta", "workspace_raw", p41 / "analysis/sv/manta_gw" / sid / "workspace/svHyGen/diploidSV.0000.vcf"),
        ("DELLY", "discover_raw", p41 / "analysis/sv/delly/discover" / f"{sid}.delly.bcf"),
        ("DELLY", "region_genotyped", p41 / "analysis/sv/delly/region_geno_20261007" / f"{sid}.geno.bcf"),
    ]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample-table", required=True, type=Path)
    ap.add_argument("--p41", required=True, type=Path)
    ap.add_argument("--out-dir", required=True, type=Path)
    ap.add_argument("--chrom", default="chr11")
    ap.add_argument("--start", type=int, default=114554000)
    ap.add_argument("--end", type=int, default=114675500)
    args = ap.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=True)
    with args.sample_table.open() as h:
        samples = list(csv.DictReader(h))
    inventory, calls = [], []
    for sample in samples:
        sid = sample["sample_id"]
        for caller, stage, path in candidate_files(args.p41, sid):
            status = {"sample_id": sid, "caller": caller, "source_stage": stage,
                      "source_path": str(path), "exists": int(path.is_file()),
                      "readable": 0, "total_records": "", "discovery_records": "", "error": ""}
            if not path.is_file():
                inventory.append(status)
                continue
            try:
                vcf = pysam.VariantFile(str(path))
                total = discovered = 0
                for rec in vcf:
                    total += 1
                    e_chrom, e_pos = endpoint(rec)
                    within = in_window(rec.chrom, rec.pos, e_chrom, e_pos, args.chrom, args.start, args.end)
                    if within:
                        discovered += 1
                        svlen = rec.info.get("SVLEN", "")
                        if isinstance(svlen, tuple):
                            svlen = ",".join(map(str, svlen))
                        calls.append({
                            "sample_id": sid, "caller": caller, "source_stage": stage,
                            "source_path": str(path), "source_readable": 1,
                            "filter": ";".join(rec.filter.keys()) or ".", "chrom": rec.chrom,
                            "pos_1based": rec.pos, "end_chrom": e_chrom, "end_1based": e_pos,
                            "svtype": rec.info.get("SVTYPE", ""), "svlen": svlen,
                            "record_id": rec.id or ".", "qual": rec.qual if rec.qual is not None else ".",
                            "alt": str(rec.alts[0]) if rec.alts else ".", "in_discovery_window": 1,
                        })
                vcf.close()
                status.update(readable=1, total_records=total, discovery_records=discovered)
            except (OSError, ValueError) as exc:
                status["error"] = f"{type(exc).__name__}: {exc}"
            inventory.append(status)
    with (args.out_dir / "caller_input_audit.tsv").open("w", newline="") as h:
        cols = ["sample_id", "caller", "source_stage", "source_path", "exists", "readable", "total_records", "discovery_records", "error"]
        w = csv.DictWriter(h, fieldnames=cols, delimiter="\t", lineterminator="\n")
        w.writeheader(); w.writerows(inventory)
    with (args.out_dir / "caller_records_in_discovery_window.tsv").open("w", newline="") as h:
        w = csv.DictWriter(h, fieldnames=FIELDS, delimiter="\t", lineterminator="\n")
        w.writeheader(); w.writerows(calls)


if __name__ == "__main__":
    main()
