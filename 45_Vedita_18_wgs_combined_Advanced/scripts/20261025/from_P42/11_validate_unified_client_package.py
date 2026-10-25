#!/usr/bin/env python3
"""Validate the client-facing P42 package without writing into that package."""
from __future__ import annotations

import argparse
import csv
import json
import re
import subprocess
from pathlib import Path


PROJECT = Path("/home/gao/projects_2026Q4/42_Vedita_18_wgs_advanced_analysis").resolve()


def args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--delivery", type=Path, required=True)
    parser.add_argument("--audit-dir", type=Path, required=True)
    parser.add_argument("--samtools", type=Path, required=True)
    parser.add_argument("--bcftools", type=Path, required=True)
    return parser.parse_args()


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def run(command: list[str]) -> None:
    subprocess.run(command, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)


def capture(command: list[str]) -> str:
    return subprocess.run(command, check=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True).stdout


def main() -> None:
    ns = args()
    delivery = ns.delivery.resolve()
    audit_dir = ns.audit_dir.resolve()
    failures: list[str] = []
    checks: dict[str, object] = {}
    if delivery.parent != PROJECT or delivery.name != "custom_research_report_20261024":
        failures.append("unexpected delivery location")
    if audit_dir.parent != PROJECT / "analysis":
        failures.append("audit directory must be a direct child of analysis/")
    audit_dir.mkdir(parents=True, exist_ok=True)
    report_files = list(delivery.glob("*.md"))
    checks["root_markdown_files"] = [path.name for path in report_files]
    if len(report_files) != 1:
        failures.append("exactly one root report is required")
    report = report_files[0] if report_files else None
    if report:
        text = report.read_text()
        required = ["## 1.", "## 2.", "## 3.", "## 4.", "## 5.", "## 6.", "## 7.", "## 8.", "## 9.", "P41 Section 7.3", "2-base-pair (bp) `AT`", "0-bp", "0/20 for Manta", "0/20 for DELLY", "907.44 GB", "not sent to the client"]
        missing = [item for item in required if item not in text]
        if missing:
            failures.append("report missing required content: " + "; ".join(missing))
        if re.search(r"long[ -]?read", text, flags=re.I):
            failures.append("report contains prohibited long-read recommendation wording")
    forbidden_dirs = [name for name in ("scripts", "logs", "analysis", "superseded_for_comparison") if (delivery / name).exists()]
    if forbidden_dirs:
        failures.append("internal directories present: " + ", ".join(forbidden_dirs))
    leaks: list[str] = []
    for suffix in ("*.md", "*.json", "*.tsv", "*.fasta"):
        for path in delivery.rglob(suffix):
            if "/mnt/" in path.read_text(errors="replace") or "/home/" in path.read_text(errors="replace"):
                leaks.append(str(path.relative_to(delivery)))
    checks["text_path_leaks"] = leaks
    if leaks:
        failures.append("internal absolute paths in client text")
    pngs = list((delivery / "figures").glob("*.png"))
    pdfs = list((delivery / "figures").glob("*.pdf"))
    checks["figures"] = {"png": len(pngs), "pdf": len(pdfs)}
    if len(pngs) != 14 or len(pdfs) != 14 or any(path.stat().st_size < 1024 for path in pngs + pdfs):
        failures.append("expected 14 non-empty PNG and 14 non-empty PDF figures")
    inventory_path = delivery / "results" / "data_inventory.tsv"
    inventory = rows(inventory_path) if inventory_path.is_file() else []
    expected = {"BAM": 18, "Manta_filtered_SV": 18, "Manta_candidate_SV": 18, "Manta_candidate_small_indel": 18, "DELLY_discovery_VCF": 18, "DELLY_region_genotyped_VCF": 18}
    observed = {key: sum(row["category"] == key for row in inventory) for key in expected}
    checks["inventory_counts"] = observed
    if observed != expected:
        failures.append("data inventory category counts differ from expected")
    snv = [row for row in inventory if row["category"] == "project_wide_SNV_indel_VCF"]
    if len(snv) != 1 or not snv[0]["availability"].startswith("unavailable"):
        failures.append("SNV/indel VCF availability declaration is missing")
    bam_total = 0
    vcf_paths: list[Path] = []
    for row in inventory:
        if row["availability"] != "included":
            continue
        path = delivery / row["client_relative_path"]
        index = delivery / row["index_client_relative_path"]
        if not path.is_file() or not index.is_file():
            failures.append("missing included data or index: " + row["client_relative_path"])
            continue
        if row["category"] == "BAM":
            bam_total += path.stat().st_size
            try:
                run([str(ns.samtools), "quickcheck", "-v", str(path)])
            except subprocess.CalledProcessError:
                failures.append("samtools quickcheck failed: " + row["client_relative_path"])
                continue
            bam_header = capture([str(ns.samtools), "view", "--no-PG", "-H", str(path)])
            if any(marker in bam_header for marker in ("/mnt/", "/home/", "/Work_bio/", "file:///")):
                failures.append("internal path remains in BAM header: " + row["client_relative_path"])
        else:
            vcf_paths.append(path)
    checks["bam_bytes"] = bam_total
    for path in vcf_paths:
        try:
            run([str(ns.bcftools), "view", "-h", str(path)])
        except subprocess.CalledProcessError:
            failures.append("bcftools header read failed: " + str(path.relative_to(delivery)))
            continue
        header = capture([str(ns.bcftools), "view", "--no-version", "-h", str(path)])
        if any(marker in header for marker in ("/mnt/", "/home/", "/Work_bio/", "file:///")):
            failures.append("internal path remains in VCF header: " + str(path.relative_to(delivery)))
    conversion = delivery / "methods" / "delly_conversion_audit.tsv"
    conversion_rows = rows(conversion) if conversion.is_file() else []
    checks["delly_conversion_rows"] = len(conversion_rows)
    if len(conversion_rows) != 36 or any(row.get("status") != "converted_and_indexed" for row in conversion_rows):
        failures.append("DELLY conversion audit is incomplete")
    bam_audit = delivery / "methods" / "bam_header_sanitization_audit.tsv"
    bam_rows = rows(bam_audit) if bam_audit.is_file() else []
    checks["bam_header_sanitization_rows"] = len(bam_rows)
    source_bam_bytes = sum(int(row.get("source_bytes", "0")) for row in bam_rows)
    checks["source_bam_bytes"] = source_bam_bytes
    if len(bam_rows) != 18 or source_bam_bytes != 907440957244 or any(row.get("source_alignment_records") != row.get("client_alignment_records") for row in bam_rows):
        failures.append("BAM header-sanitization audit is incomplete or changed alignment-record counts")
    header_audit = delivery / "methods" / "vcf_header_sanitization_audit.tsv"
    header_rows = rows(header_audit) if header_audit.is_file() else []
    checks["vcf_header_sanitization_rows"] = len(header_rows)
    if len(header_rows) != 90 or any(row.get("source_variant_body_sha256") != row.get("client_variant_body_sha256") for row in header_rows):
        failures.append("VCF header-sanitization audit is incomplete or changed variant records")
    required_results = ["c2_junction_resolution_summary.tsv", "c2_read_end_junction_resolution.tsv", "candidate_manta_delly_concordance.tsv", "clone_specific_comparison_with_criteria.tsv", "read_rescue_candidate_impact.tsv", "assay_design_selection.tsv"]
    missing_results = [name for name in required_results if not (delivery / "results" / name).is_file()]
    if missing_results:
        failures.append("missing required result tables: " + ", ".join(missing_results))
    caller = rows(delivery / "results" / "candidate_manta_delly_concordance.tsv") if (delivery / "results" / "candidate_manta_delly_concordance.tsv").is_file() else []
    checks["caller_table_rows"] = len(caller)
    if not caller or not any("0" in " ".join(row.values()) for row in caller):
        failures.append("caller table lacks expected zero-concordance records")
    payload = {"delivery": str(delivery), "pass": not failures, "failures": failures, "checks": checks}
    (audit_dir / "client_package_self_audit.json").write_text(json.dumps(payload, indent=2) + "\n")
    if failures:
        raise SystemExit("Validation failed: " + " | ".join(failures))
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
