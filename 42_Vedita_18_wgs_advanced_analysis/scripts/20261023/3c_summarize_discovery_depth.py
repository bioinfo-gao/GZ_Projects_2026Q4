#!/usr/bin/env python3
"""Summarize audited per-base P42 discovery depth into fixed 1-kb bins.

Input tracks were generated with identical FLAG handling; only MAPQ threshold
differs.  Zero-coverage positions are retained, so mean depth and zero fraction
remain interpretable in repeat/low-mappability intervals.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import statistics
from pathlib import Path

TRACKS = {
    "standard_primary_mapq20": "standard_primary_mapq20.depth.tsv.gz",
    "inclusive_mapq0": "inclusive_mapq0.depth.tsv.gz",
}


def summarize(path: Path, sample_id: str, track: str, bin_size: int, start_1based: int):
    bins: dict[int, list[int]] = {}
    chrom_seen: set[str] = set()
    with gzip.open(path, "rt") as handle:
        for line in handle:
            chrom, position, depth = line.rstrip("\n").split("\t")
            pos = int(position)
            chrom_seen.add(chrom)
            bin_start = start_1based + ((pos - start_1based) // bin_size) * bin_size
            bins.setdefault(bin_start, []).append(int(depth))
    if chrom_seen != {"chr11"}:
        raise ValueError(f"Unexpected chromosome(s) in {path}: {sorted(chrom_seen)}")
    rows = []
    for bin_start, values in sorted(bins.items()):
        rows.append({
            "sample_id": sample_id, "track": track, "chrom": "chr11",
            "start_1based": bin_start, "end_1based": bin_start + len(values) - 1,
            "n_bases": len(values), "mean_depth": f"{statistics.fmean(values):.6f}",
            "median_depth": f"{statistics.median(values):.6f}",
            "min_depth": min(values), "max_depth": max(values),
            "zero_depth_fraction": f"{sum(v == 0 for v in values) / len(values):.6f}",
        })
    return rows


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--input-dir", required=True, type=Path)
    p.add_argument("--out", required=True, type=Path)
    p.add_argument("--start", type=int, default=114554000)
    p.add_argument("--bin-size", type=int, default=1000)
    args = p.parse_args()
    output_rows = []
    sample_dirs = sorted(path for path in args.input_dir.iterdir() if path.is_dir())
    if len(sample_dirs) != 18:
        raise SystemExit(f"Expected exactly 18 sample directories, found {len(sample_dirs)}")
    for sample_dir in sample_dirs:
        for track, filename in TRACKS.items():
            path = sample_dir / filename
            if not path.is_file():
                raise SystemExit(f"Missing required depth file: {path}")
            output_rows.extend(summarize(path, sample_dir.name, track, args.bin_size, args.start))
    args.out.parent.mkdir(parents=True, exist_ok=True)
    fields = ["sample_id", "track", "chrom", "start_1based", "end_1based", "n_bases",
              "mean_depth", "median_depth", "min_depth", "max_depth", "zero_depth_fraction"]
    with args.out.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows(output_rows)


if __name__ == "__main__":
    main()
