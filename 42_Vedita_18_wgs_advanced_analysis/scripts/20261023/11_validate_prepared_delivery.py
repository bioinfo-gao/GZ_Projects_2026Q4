#!/usr/bin/env python3
"""Fail-closed self-audit of P42's prepared, unsent client delivery."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
from collections import Counter
from pathlib import Path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def tsv_rows(path: Path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--project", type=Path, required=True)
    parser.add_argument("--delivery", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    args = parser.parse_args()
    args.out_dir.mkdir(parents=True, exist_ok=False)
    failures = []
    project = args.project
    delivery = args.delivery
    report = delivery / "P42_Vedita_18_WGS_1023.md"
    if delivery.name != "custom_research_report_20261023": failures.append("delivery folder date is not 20261023")
    if not report.is_file(): failures.append("expected single report file is missing")
    reports = list(delivery.glob("*.md"))
    if reports != [report]: failures.append(f"expected exactly one top-level report, observed {[p.name for p in reports]}")
    text = report.read_text() if report.is_file() else ""
    if "**Report Date:** 2026-10-23" not in text: failures.append("report body date mismatch")
    required_sections = [
        "## 1. Objectives", "## 2. Core Conclusions", "## 3. Key Findings", "## 4. Sample Information",
        "## 5. Analysis Rationale and Decision Criteria", "## 6. Methods", "## 7. Results",
        "## 8. Conclusions, Limitations and Next Steps", "### 8.2 Practical reading with evidence sources", "## 9. Deliverable Files",
    ]
    for section in required_sections:
        if section not in text: failures.append(f"missing required report section: {section}")
    for prohibited in (r"the same (answers|ground|findings)", r"repeats the", r"restates the"):
        if re.search(prohibited, text, flags=re.IGNORECASE): failures.append(f"simple-repetition wording found: {prohibited}")
    for required in ("Zhen Gao, PhD", "Principal Bioinformatics Scientist", "Athenomics", "Tissue / Cell:** not specified by client"):
        if required not in text: failures.append(f"missing required report header/signature content: {required}")
    expected_figures = sorted((delivery / "figures").glob("P42_Fig*.png"))
    expected_pdfs = sorted((delivery / "figures").glob("P42_Fig*.pdf"))
    if len(expected_figures) != 12 or len(expected_pdfs) != 12:
        failures.append(f"expected 12 PNG and 12 PDF figures, observed {len(expected_figures)} PNG and {len(expected_pdfs)} PDF")
    manifest_path = delivery / "figures/figure_manifest.json"
    manifest = json.loads(manifest_path.read_text()) if manifest_path.is_file() else {}
    if manifest.get("figure_count") != 12: failures.append("figure manifest count mismatch")
    if manifest.get("case_only_candidate_count") != 4: failures.append("case-only candidate count mismatch")
    if manifest.get("high_confidence_dox_candidate_count") != 0: failures.append("high-confidence candidate count mismatch")

    out_root = project / "analysis/20261023_fullcohort_repeat_resolve"
    copies = {
        "results/fixed_candidate_summary.tsv": out_root / "09_fixed_candidate_regenotyping_v1/fixed_candidate_summary.tsv",
        "results/fixed_candidate_dld1_matrix.tsv": out_root / "09_fixed_candidate_regenotyping_v1/fixed_candidate_dld1_matrix.tsv",
        "results/raw_evidence_match_audit.tsv.gz": out_root / "09_fixed_candidate_regenotyping_v1/raw_evidence_match_audit.tsv.gz",
        "results/contig_placement_summary.tsv": out_root / "05_local_assembly/contig_placements_v1/contig_placement_summary.tsv",
        "methods/analysis_audit.json": out_root / "09_fixed_candidate_regenotyping_v1/audit.json",
    }
    copy_hashes = {}
    for relative, source in copies.items():
        target = delivery / relative
        if not source.is_file() or not target.is_file():
            failures.append(f"missing source or delivery copy for {relative}")
            continue
        copy_hashes[relative] = {"source": sha256(source), "delivery": sha256(target)}
        if copy_hashes[relative]["source"] != copy_hashes[relative]["delivery"]:
            failures.append(f"delivery copy differs from source: {relative}")

    summary = tsv_rows(delivery / "results/fixed_candidate_summary.tsv")
    matrix = tsv_rows(delivery / "results/fixed_candidate_dld1_matrix.tsv")
    classifications = Counter(row["classification"] for row in summary)
    if len(summary) != 282: failures.append(f"expected 282 fixed candidates, observed {len(summary)}")
    if len(matrix) != 1410: failures.append(f"expected 1410 candidate-sample matrix rows, observed {len(matrix)}")
    expected_classifications = {
        "anchor_context_without_matching_dual_endpoint": 62,
        "case_only_raw_candidate_requires_sensitivity_and_assembly": 4,
        "control_only_raw_dual_endpoint_evidence": 21,
        "not_evaluable_as_two_endpoint_bnd": 188,
        "shared_raw_dual_endpoint_evidence": 7,
    }
    if dict(classifications) != expected_classifications:
        failures.append(f"classification totals mismatch: {dict(classifications)}")
    contigs = tsv_rows(delivery / "results/contig_placement_summary.tsv")
    if len(contigs) != 90: failures.append(f"expected 90 contigs, observed {len(contigs)}")
    if sum(int(row["distinct_loci"] or 0) > 1 for row in contigs) != 32:
        failures.append("alternative-placement contig total mismatch")
    claims = ["282 caller breakend events", "94 had a uniquely parsed two-endpoint geometry", "four case-only candidate patterns", "32 of 90"]
    for claim in claims:
        if claim not in text: failures.append(f"numeric report claim missing: {claim}")
    audit = {
        "report": str(report), "report_date": "2026-10-23", "delivery_folder": delivery.name,
        "figure_png_count": len(expected_figures), "figure_pdf_count": len(expected_pdfs),
        "fixed_candidate_count": len(summary), "candidate_sample_rows": len(matrix),
        "classification_counts": dict(sorted(classifications.items())),
        "contig_count": len(contigs), "alternative_placement_contigs": sum(int(row["distinct_loci"] or 0) > 1 for row in contigs),
        "copy_sha256": copy_hashes, "status": "pass" if not failures else "fail", "failures": failures,
    }
    (args.out_dir / "delivery_self_audit.json").write_text(json.dumps(audit, indent=2, sort_keys=True) + "\n")
    if failures:
        raise SystemExit("P42 prepared-delivery self-audit failed: " + "; ".join(failures))
    (args.out_dir / ".complete").write_text("completed\n")
    print(json.dumps(audit, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
