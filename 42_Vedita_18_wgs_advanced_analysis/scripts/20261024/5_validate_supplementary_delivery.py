#!/usr/bin/env python3
"""Self-audit the P42 2026-10-24 incremental delivery without mutating it."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
from datetime import datetime
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--delivery", type=Path, required=True)
    parser.add_argument("--audit-dir", type=Path, required=True)
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def tsv_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def check(condition: bool, message: str, checks: list[dict[str, object]]) -> None:
    checks.append({"check": message, "passed": condition})
    if not condition:
        raise AssertionError(message)


def main() -> None:
    args = parse_args()
    out_dir = args.out_dir.resolve()
    delivery = args.delivery.resolve()
    audit = args.audit_dir.resolve()
    if audit.parent != out_dir:
        raise ValueError(f"Audit directory must be a direct child of out-dir: {audit}")
    if audit.exists():
        raise FileExistsError(f"Refusing to overwrite prior audit directory: {audit}")
    audit.mkdir()
    checks: list[dict[str, object]] = []

    try:
        report = delivery / "P42_P41_7_3_Assay_Followup_1024.md"
        check(report.is_file(), "English supplementary report exists", checks)
        report_text = report.read_text()
        for section in range(1, 10):
            check(bool(re.search(rf"^## {section}\. ", report_text, flags=re.MULTILINE)), f"Report contains Section {section}", checks)
        check("P41" in report_text and "Section 7.3" in report_text, "Report explicitly states P41 Section 7.3 continuity", checks)
        check("No PCR, new sequencing, or raw FISH" in report_text, "Report states external-data limitation", checks)
        check("long-read" not in report_text.lower() and "long read" not in report_text.lower(), "Report contains no prohibited long-read recommendation", checks)
        check("Prepared for scientific review; not sent to the client." in report_text, "Report status is review-only", checks)

        figures = delivery / "figures"
        pngs = sorted(figures.glob("*.png"))
        pdfs = sorted(figures.glob("*.pdf"))
        check(len(pngs) == 1, "Exactly one PNG figure supplied", checks)
        check(len(pdfs) == 1, "Exactly one PDF figure supplied", checks)
        check(pngs[0].stem == pdfs[0].stem, "PNG/PDF figure basenames match", checks)
        check(pngs[0].stat().st_size > 10_000 and pdfs[0].stat().st_size > 10_000, "Figure files are nontrivial", checks)

        selection_source = out_dir / "assay_design_selection.tsv"
        selection_delivery = delivery / "results" / "assay_design_selection.tsv"
        check(selection_delivery.is_file(), "Assay selection table supplied", checks)
        selection = tsv_rows(selection_delivery)
        check({row["candidate_label"] for row in selection} == {"C1", "C2", "C3", "C4"}, "Selection table contains C1–C4 exactly", checks)
        check(sha256(selection_source) == sha256(selection_delivery), "Delivered selection table matches analysis source", checks)

        for name in ["primer_specificity.tsv", "candidate_template_audit.tsv", "raw_support_alignment_audit.tsv"]:
            source = out_dir / name
            copied = delivery / "results" / name
            check(copied.is_file(), f"{name} supplied", checks)
            check(sha256(source) == sha256(copied), f"{name} hash matches analysis source", checks)
        template_source = out_dir / "virtual_junction_templates.fasta"
        template_copy = delivery / "methods" / "virtual_junction_templates.fasta"
        check(template_copy.is_file() and sha256(template_source) == sha256(template_copy), "Virtual-template FASTA hash matches analysis source", checks)

        external = tsv_rows(delivery / "results" / "external_input_status.tsv")
        check(len(external) == 3, "External-input status table has three required evidence categories", checks)
        check(all(row["status"] == "not_provided" for row in external), "External-input status does not claim absent data were analysed", checks)

        forbidden_extensions = {".bam", ".cram", ".fastq", ".fq", ".vcf", ".bcf"}
        disallowed = [str(path.relative_to(delivery)) for path in delivery.rglob("*") if path.is_file() and path.suffix.lower() in forbidden_extensions]
        check(not disallowed, "Supplement contains no duplicated large sequencing files", checks)
        check(not (delivery / "P42_Vedita_18_WGS_1023.md").exists(), "Supplement does not duplicate the prior P42 report", checks)

        success = True
        error = None
    except Exception as exc:
        success = False
        error = str(exc)

    payload = {
        "timestamp": datetime.now().astimezone().isoformat(timespec="seconds"),
        "delivery": str(delivery),
        "passed": success,
        "checks": checks,
        "error": error,
    }
    (audit / "delivery_self_audit.json").write_text(json.dumps(payload, indent=2) + "\n")
    if not success:
        raise SystemExit(error)


if __name__ == "__main__":
    main()
