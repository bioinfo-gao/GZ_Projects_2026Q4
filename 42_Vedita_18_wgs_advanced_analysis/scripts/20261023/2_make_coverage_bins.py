#!/usr/bin/env python3
"""Create reproducible, non-overlapping coverage bins from the exact reference fai."""
from __future__ import annotations
import argparse
from pathlib import Path

CANONICAL = {f"chr{i}" for i in range(1, 23)} | {"chrX", "chrY"}

def write_bins(path: Path, contigs, width: int):
    with path.open("w") as out:
        for chrom, length in contigs:
            for start in range(0, length, width):
                out.write(f"{chrom}\t{start}\t{min(start + width, length)}\n")

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--fai", required=True, type=Path)
    p.add_argument("--out-dir", required=True, type=Path)
    args = p.parse_args(); args.out_dir.mkdir(parents=True, exist_ok=True)
    lengths = {}
    for line in args.fai.read_text().splitlines():
        fields = line.split("\t")
        lengths[fields[0]] = int(fields[1])
    canonical = [(c, lengths[c]) for c in lengths if c in CANONICAL]
    if len(canonical) != 24 or "chr11" not in lengths:
        raise SystemExit("Expected canonical chr1-22, chrX, chrY and chr11 in reference fai")
    write_bins(args.out_dir / "canonical_100kb.bed", canonical, 100_000)
    write_bins(args.out_dir / "chr11_10kb.bed", [("chr11", lengths["chr11"])], 10_000)
    with (args.out_dir / "discovery_1kb.bed").open("w") as out:
        for start in range(114_554_000 - 1, 114_675_500, 1000):
            out.write(f"chr11\t{start}\t{min(start + 1000, 114_675_500)}\n")

if __name__ == "__main__":
    main()
