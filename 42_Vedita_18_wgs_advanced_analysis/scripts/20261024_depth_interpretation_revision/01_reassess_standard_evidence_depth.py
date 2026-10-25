#!/usr/bin/env python3
"""Recalculate the analysis-specific autosomal depth context without recasting it as service QC.

This revision deliberately separates the WGS data-yield context (approximately 25-30x for most
samples, as supplied by the project owner) from the restricted denominator used in the C1-C4
comparison: primary, non-duplicate, non-QC-fail, non-supplementary alignments with MAPQ >=20.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import json
import math
from pathlib import Path


AUTOSOMES = {f"chr{number}" for number in range(1, 23)}
DLD1_ROLES = {
    "S01_DLD1_pop": "DLD1 parental population",
    "S02_DLD1_c1_ctrl": "DLD1 control clone",
    "S03_DLD1_c2_ctrl": "DLD1 control clone",
    "S04_DLD1_c3_dox": "DLD1 +doxycycline clone",
    "S05_DLD1_c4_dox": "DLD1 +doxycycline clone",
}
DEPTH_DEFINITION = (
    "Base-weighted chr1-chr22 depth from the existing standard evidence track: "
    "primary, non-duplicate, non-QC-fail, non-supplementary alignments with MAPQ >=20. "
    "This is an analysis-specific denominator, not the WGS service-yield metric."
)


def arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--coverage-root", required=True, type=Path)
    parser.add_argument("--out-dir", required=True, type=Path)
    return parser.parse_args()


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def autosomal_mean(regions_path: Path) -> tuple[float, int]:
    bases = 0
    depth_sum = 0.0
    with gzip.open(regions_path, "rt") as handle:
        for line in handle:
            chrom, start, end, depth = line.rstrip("\n").split("\t")
            if chrom not in AUTOSOMES:
                continue
            width = int(end) - int(start)
            bases += width
            depth_sum += width * float(depth)
    if not bases:
        raise ValueError(f"No autosomal bins in {regions_path}")
    return depth_sum / bases, bases


def probability_at_least(depth: float, fraction: float, count: int) -> float:
    n = max(0, round(depth))
    return 1.0 - sum(
        math.comb(n, k) * fraction**k * (1.0 - fraction) ** (n - k)
        for k in range(count)
    )


def main() -> None:
    ns = arguments()
    coverage_root = ns.coverage_root.resolve()
    out_dir = ns.out_dir.resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, object]] = []
    for sample_dir in sorted(path for path in coverage_root.iterdir() if path.is_dir()):
        files = list(sample_dir.glob("*.regions.bed.gz"))
        if len(files) != 1:
            raise ValueError(f"Expected exactly one regions file in {sample_dir}, found {len(files)}")
        depth, bases = autosomal_mean(files[0])
        sample_id = sample_dir.name
        rows.append(
            {
                "sample_id": sample_id,
                "analysis_group": DLD1_ROLES.get(sample_id, "non-DLD1 cohort background"),
                "autosomal_standard_evidence_mean_depth": f"{depth:.6f}",
                "autosomal_bases": bases,
                "depth_definition": DEPTH_DEFINITION,
                "data_yield_context": "Approximately 25-30x WGS data amount for most samples; separate from this restricted analysis denominator.",
                "interpretation": "Use only to bound an absence-style interpretation in the standard evidence comparison; do not use as a service-quality or pass/fail metric.",
            }
        )
    if len(rows) != 18:
        raise ValueError(f"Expected 18 samples, found {len(rows)}")
    base_counts = {int(row["autosomal_bases"]) for row in rows}
    if len(base_counts) != 1:
        raise ValueError(f"Autosomal base count differs among samples: {sorted(base_counts)}")
    fields = list(rows[0])
    write_tsv(out_dir / "standard_evidence_autosomal_depth.tsv", rows, fields)

    dld1 = [row for row in rows if row["sample_id"] in DLD1_ROLES]
    depths = [float(row["autosomal_standard_evidence_mean_depth"]) for row in rows]
    dld1_depths = [float(row["autosomal_standard_evidence_mean_depth"]) for row in dld1]
    summary = [
        {
            "scope": "All 18 supplied WGS samples",
            "sample_count": len(rows),
            "minimum_standard_evidence_depth": f"{min(depths):.6f}",
            "maximum_standard_evidence_depth": f"{max(depths):.6f}",
            "depth_definition": DEPTH_DEFINITION,
            "interpretation": "Describes the restricted denominator used in this analysis, not total WGS data yield or service conformance.",
        },
        {
            "scope": "DLD1 S01-S05 clone comparison",
            "sample_count": len(dld1),
            "minimum_standard_evidence_depth": f"{min(dld1_depths):.6f}",
            "maximum_standard_evidence_depth": f"{max(dld1_depths):.6f}",
            "depth_definition": DEPTH_DEFINITION,
            "interpretation": "Comparator zeros are reported as not observed within this evidence context; no retrospective numeric depth pass/fail gate is applied.",
        },
    ]
    write_tsv(out_dir / "standard_evidence_depth_summary.tsv", summary, list(summary[0]))

    nominal_rows: list[dict[str, object]] = []
    for row in dld1:
        depth = float(row["autosomal_standard_evidence_mean_depth"])
        nominal_rows.append(
            {
                "sample_id": row["sample_id"],
                "standard_evidence_mean_depth": f"{depth:.6f}",
                "nominal_event_fraction": "0.05",
                "probability_at_least_one_supporting_fragment": f"{probability_at_least(depth, 0.05, 1):.6f}",
                "probability_at_least_three_supporting_fragments": f"{probability_at_least(depth, 0.05, 3):.6f}",
                "model_boundary": "Illustrative independent-binomial calibration only. It is not an allele-fraction estimate, an assay validation, a repeat-locus sensitivity claim, or a decision threshold.",
            }
        )
    write_tsv(out_dir / "nominal_5pct_detection_context.tsv", nominal_rows, list(nominal_rows[0]))

    audit = {
        "sample_count": len(rows),
        "autosomal_base_count_per_sample": next(iter(base_counts)),
        "dld1_sample_count": len(dld1),
        "all_sample_depth_range": [round(min(depths), 6), round(max(depths), 6)],
        "dld1_depth_range": [round(min(dld1_depths), 6), round(max(dld1_depths), 6)],
        "inapplicable_rules_removed": ["No retrospective numeric depth gate", "No direct support/depth allele-fraction estimate", "No absence conclusion from a comparator zero"],
        "pass": True,
    }
    (out_dir / "depth_interpretation_audit.json").write_text(json.dumps(audit, indent=2) + "\n")
    print(json.dumps(audit, indent=2))


if __name__ == "__main__":
    main()
