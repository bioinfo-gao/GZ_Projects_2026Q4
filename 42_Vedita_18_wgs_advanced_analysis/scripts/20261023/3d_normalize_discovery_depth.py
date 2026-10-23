#!/usr/bin/env python3
"""Normalize local 1-kb depth by each sample's canonical-genome mean depth."""
from __future__ import annotations

import argparse
import csv
import gzip
from pathlib import Path


def genome_means(coverage_root: Path) -> dict[str, float]:
    means = {}
    for sample_dir in sorted(path for path in coverage_root.iterdir() if path.is_dir()):
        files = list(sample_dir.glob("*.regions.bed.gz"))
        if len(files) != 1:
            raise ValueError(f"Expected one genome regions file in {sample_dir}, found {len(files)}")
        total_bases = 0; total_depth = 0.0
        with gzip.open(files[0], "rt") as handle:
            for line in handle:
                chrom, start, end, depth = line.rstrip("\n").split("\t")
                width = int(end) - int(start)
                total_bases += width; total_depth += width * float(depth)
        if not total_bases:
            raise ValueError(f"No canonical bins in {files[0]}")
        means[sample_dir.name] = total_depth / total_bases
    if len(means) != 18:
        raise ValueError(f"Expected 18 genome means, found {len(means)}")
    return means


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--summary", required=True, type=Path)
    p.add_argument("--genome-coverage-dir", required=True, type=Path)
    p.add_argument("--out", required=True, type=Path)
    args = p.parse_args()
    means = genome_means(args.genome_coverage_dir)
    with args.summary.open(newline="") as source:
        rows = list(csv.DictReader(source, delimiter="\t"))
        fields = list(rows[0]) + ["canonical_genome_mean_depth", "depth_over_genome_mean"]
    for row in rows:
        mean = means[row["sample_id"]]
        row["canonical_genome_mean_depth"] = f"{mean:.6f}"
        row["depth_over_genome_mean"] = f"{float(row['mean_depth']) / mean:.6f}"
    with args.out.open("w", newline="") as target:
        writer = csv.DictWriter(target, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(rows)
    with args.out.with_name("canonical_genome_mean_depth.tsv").open("w", newline="") as target:
        writer = csv.writer(target, delimiter="\t", lineterminator="\n")
        writer.writerow(["sample_id", "canonical_genome_mean_depth"])
        for sample_id, mean in sorted(means.items()): writer.writerow([sample_id, f"{mean:.6f}"])


if __name__ == "__main__":
    main()
