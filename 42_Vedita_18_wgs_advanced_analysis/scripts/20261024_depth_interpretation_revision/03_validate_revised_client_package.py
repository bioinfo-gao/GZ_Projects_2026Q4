#!/usr/bin/env python3
"""Validate the revised client package without modifying it.

The unchanged BAM/VCF/figure payload is hard-linked from the archived package that already passed
the full record-level self-audit. This validator proves the link identity, rechecks readability and
header privacy, and separately checks the revised report and result-table interpretation.
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import re
import subprocess
from pathlib import Path


EXPECTED = {
    "BAM": 18,
    "Manta_filtered_SV": 18,
    "Manta_candidate_SV": 18,
    "Manta_candidate_small_indel": 18,
    "DELLY_discovery_VCF": 18,
    "DELLY_region_genotyped_VCF": 18,
}
PATH_MARKERS = ("/mnt/", "/home/", "/Work_bio/", "file:///")
FORBIDDEN_INTERPRETATION = ("60x", "60×", "service nonconformance", "below service target")


def arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--source-package", required=True, type=Path)
    parser.add_argument("--delivery", required=True, type=Path)
    parser.add_argument("--audit-dir", required=True, type=Path)
    parser.add_argument("--samtools", required=True, type=Path)
    parser.add_argument("--bcftools", required=True, type=Path)
    return parser.parse_args()


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def capture(command: list[str]) -> str:
    return subprocess.run(command, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True).stdout


def run(command: list[str]) -> None:
    subprocess.run(command, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)


def main() -> None:
    ns = arguments()
    source = ns.source_package.resolve()
    delivery = ns.delivery.resolve()
    audit_dir = ns.audit_dir.resolve()
    failures: list[str] = []
    checks: dict[str, object] = {}
    audit_dir.mkdir(parents=True, exist_ok=True)
    report_files = list(delivery.glob("*.md"))
    if len(report_files) != 1:
        failures.append("exactly one root report is required")
        report = None
    else:
        report = report_files[0]
        text = report.read_text()
        required = [
            "## 1.", "## 2.", "## 3.", "## 4.", "## 5.", "## 6.", "## 7.", "## 8.", "## 9.",
            "P41 Section 7.3", "0/20 for Manta", "0/20 for DELLY", "2-base-pair (bp) `AT`", "0-bp",
            "approximately 25-30x", "not sent to the client", "not as the WGS service-yield metric",
        ]
        missing = [entry for entry in required if entry not in text]
        if missing:
            failures.append("report missing required content: " + "; ".join(missing))
    checks["root_markdown_files"] = [path.name for path in report_files]
    forbidden_dirs = [name for name in ("scripts", "logs", "analysis", "superseded_for_comparison") if (delivery / name).exists()]
    if forbidden_dirs:
        failures.append("internal directories found in client package: " + ", ".join(forbidden_dirs))

    text_leaks: list[str] = []
    interpretation_leaks: list[str] = []
    for suffix in ("*.md", "*.json", "*.tsv", "*.fasta", "*.txt"):
        for path in delivery.rglob(suffix):
            text = path.read_text(errors="replace")
            if any(marker in text for marker in PATH_MARKERS):
                text_leaks.append(str(path.relative_to(delivery)))
            if any(marker.lower() in text.lower() for marker in FORBIDDEN_INTERPRETATION):
                interpretation_leaks.append(str(path.relative_to(delivery)))
    checks["text_path_leaks"] = text_leaks
    checks["inapplicable_depth_or_service_wording"] = interpretation_leaks
    if text_leaks:
        failures.append("internal path leak in client text")
    if interpretation_leaks:
        failures.append("inapplicable depth-gate or service-quality wording remains")

    manifest_path = delivery / "figures" / "figure_manifest.json"
    manifest = json.loads(manifest_path.read_text()) if manifest_path.is_file() else {}
    pngs = list((delivery / "figures").glob("*.png"))
    pdfs = list((delivery / "figures").glob("*.pdf"))
    checks["figures"] = {"manifest_count": manifest.get("figure_count"), "png": len(pngs), "pdf": len(pdfs)}
    if manifest.get("figure_count") != 14 or len(pngs) != 14 or len(pdfs) != 14 or any(path.stat().st_size < 1024 for path in pngs + pdfs):
        failures.append("14 non-empty PNG and 14 non-empty PDF figures are required")
    fig12_pdf = delivery / "figures" / "P42_Fig12_evidence_boundary_summary.pdf"
    try:
        fig12_text = capture(["pdftotext", str(fig12_pdf), "-"])
        required_fig12_text = (
            "Clone/doxycycline specificity is not established from the current WGS comparison alone",
            "Comparator zeros are reported as not observed",
            "not a service-yield metric",
        )
        missing_fig12_text = [item for item in required_fig12_text if item not in fig12_text]
        prohibited_fig12_text = [item for item in ("60", "NOT MET", "below the") if item in fig12_text]
        checks["figure12_pdf_text"] = {"missing_required": missing_fig12_text, "prohibited": prohibited_fig12_text}
        if missing_fig12_text or prohibited_fig12_text:
            failures.append("Figure 12 retains or omits critical depth-interpretation wording")
    except (FileNotFoundError, subprocess.CalledProcessError):
        failures.append("Figure 12 PDF text could not be checked")

    inventory = rows(delivery / "results" / "data_inventory.tsv")
    observed = {category: sum(row["category"] == category for row in inventory) for category in EXPECTED}
    checks["inventory_counts"] = observed
    if observed != EXPECTED:
        failures.append("inventory counts differ from expected")
    inode_mismatches: list[str] = []
    bam_checked = 0
    vcf_checked = 0
    for row in inventory:
        if row["availability"] != "included":
            continue
        relative = Path(row["client_relative_path"])
        client_path = delivery / relative
        source_path = source / relative
        client_index = delivery / row["index_client_relative_path"]
        source_index = source / row["index_client_relative_path"]
        if not client_path.is_file() or not source_path.is_file() or not client_index.is_file() or not source_index.is_file():
            failures.append(f"missing data or index: {relative}")
            continue
        if os.stat(client_path).st_ino != os.stat(source_path).st_ino or os.stat(client_index).st_ino != os.stat(source_index).st_ino:
            inode_mismatches.append(str(relative))
        try:
            if row["category"] == "BAM":
                run([str(ns.samtools), "quickcheck", "-v", str(client_path)])
                header = capture([str(ns.samtools), "view", "--no-PG", "-H", str(client_path)])
                bam_checked += 1
            else:
                header = capture([str(ns.bcftools), "view", "--no-version", "-h", str(client_path)])
                vcf_checked += 1
            if any(marker in header for marker in PATH_MARKERS):
                failures.append(f"internal path remains in header: {relative}")
        except subprocess.CalledProcessError:
            failures.append(f"readability check failed: {relative}")
    checks["hard_link_inode_mismatches"] = inode_mismatches
    checks["bam_quickcheck_count"] = bam_checked
    checks["vcf_header_check_count"] = vcf_checked
    if inode_mismatches:
        failures.append("unchanged data payload is not identity-linked to the self-audited archive")

    clone_rows = rows(delivery / "results" / "clone_specific_comparison_with_criteria.tsv")
    depth_rows = rows(delivery / "results" / "standard_evidence_autosomal_depth.tsv")
    depth_summary = rows(delivery / "results" / "standard_evidence_depth_summary.tsv")
    checks["revised_result_rows"] = {"clone": len(clone_rows), "depth": len(depth_rows), "depth_summary": len(depth_summary)}
    if len(clone_rows) != 20 or len(depth_rows) != 18 or len(depth_summary) != 2:
        failures.append("revised results have an unexpected row count")
    if any("60" in row.get("filtering_rule", "") for row in clone_rows):
        failures.append("clone table retains an inapplicable numeric depth gate")
    if not all("not the WGS service-yield metric" in row.get("depth_definition", "") for row in depth_rows):
        failures.append("depth table does not clearly separate analysis denominator from service yield")
    client_method = rows(delivery / "methods" / "depth_interpretation_method.tsv")
    if len(client_method) != 3:
        failures.append("client depth-method table is missing or incomplete")

    payload = {"delivery": str(delivery), "source_archive": str(source), "pass": not failures, "failures": failures, "checks": checks}
    (audit_dir / "client_package_self_audit.json").write_text(json.dumps(payload, indent=2) + "\n")
    if failures:
        raise SystemExit("Validation failed: " + " | ".join(failures))
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
