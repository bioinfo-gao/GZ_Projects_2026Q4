#!/usr/bin/env python3
"""Self-audit the P42 client-question supplementary package without modifying it."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
from datetime import datetime
from pathlib import Path


REPORT_NAME = "P42_Client_Request_Completion_1024.md"
FIGURE_STEM = "P42_Fig14_Client_Request_Completion"


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
    checks.append({"check": message, "passed": bool(condition)})
    if not condition:
        raise AssertionError(message)


def main() -> None:
    args = parse_args()
    out_dir = args.out_dir.resolve()
    delivery = args.delivery.resolve()
    audit = args.audit_dir.resolve()
    if audit.parent != out_dir:
        raise ValueError(f"Audit directory must be a direct child of analysis output: {audit}")
    if audit.exists():
        raise FileExistsError(f"Refusing to overwrite a prior audit: {audit}")
    audit.mkdir()
    checks: list[dict[str, object]] = []
    try:
        report = delivery / REPORT_NAME
        check(report.is_file(), "English client-review report exists", checks)
        report_text = report.read_text()
        for section in range(1, 10):
            check(bool(re.search(rf"^## {section}\. ", report_text, flags=re.MULTILINE)),
                  f"Report contains mandatory Section {section}", checks)
        check("### 8.2 Practical reading with evidence paths" in report_text,
              "Report includes the evidence-linked Section 8.2 practical reading", checks)
        check("P41 Section 7.3" in report_text and "C2" in report_text,
              "Report explicitly states P41-to-P42 continuity", checks)
        check("Zhen Gao, PhD" in report_text and "Athenomics" in report_text,
              "Report carries required author and organization signature", checks)
        check("Analysis Platform:" in report_text and "Species:" in report_text and "Tissue / Cell:" in report_text,
              "Report header includes platform, species, and cell context", checks)
        check("2-bp `AT`" in report_text and "0-bp" in report_text,
              "Report states the C2 microhomology and gap result", checks)
        check("0/20 for Manta" in report_text and "0/20 for DELLY" in report_text and "0/20 for both callers" in report_text,
              "Report states exact caller-concordance totals", checks)
        check("907.44 GB" in report_text and "SNV/indel VCF" in report_text,
              "Report states BAM volume and SNV/indel availability boundary", checks)
        check("long-read" not in report_text.lower() and "long read" not in report_text.lower(),
              "Report contains no prohibited long-read recommendation", checks)
        check("Prepared for scientific review; not sent to the client." in report_text,
              "Report is explicitly review-only", checks)
        check(not re.search(r"the same (answers|ground|findings)|repeats the|restates the", report_text, re.I),
              "Report does not frame layered sections as simple repetition", checks)

        figures = delivery / "figures"
        png = figures / f"{FIGURE_STEM}.png"
        pdf = figures / f"{FIGURE_STEM}.pdf"
        check(png.is_file() and pdf.is_file(), "Both PNG and PDF Figure 14 files are supplied", checks)
        check(png.stat().st_size > 10_000 and pdf.stat().st_size > 10_000,
              "Figure files are nontrivial", checks)
        check(len(list(figures.glob("*.png"))) == 1 and len(list(figures.glob("*.pdf"))) == 1,
              "Supplement has exactly one PNG/PDF figure pair", checks)

        result_names = [
            "c2_junction_resolution_summary.tsv",
            "c2_read_end_junction_resolution.tsv",
            "candidate_manta_delly_concordance.tsv",
            "clone_specific_comparison_with_criteria.tsv",
            "read_rescue_candidate_impact.tsv",
            "customer_request_coverage.tsv",
            "data_delivery_availability_manifest.tsv",
            "customer_request_audit.json",
        ]
        for name in result_names:
            source = out_dir / name
            copied = delivery / "results" / name
            check(copied.is_file(), f"New result {name} is supplied", checks)
            check(sha256(source) == sha256(copied), f"{name} hash matches its analysis source", checks)

        c2 = tsv_rows(delivery / "results" / "c2_junction_resolution_summary.tsv")
        check(len(c2) == 1 and c2[0]["independent_qname_count"] == "3" and c2[0]["split_read_end_count"] == "5",
              "Delivered C2 summary has three QNAMEs and five read ends", checks)
        check(c2[0]["query_overlap_bp_distribution"] == "2bp:5" and c2[0]["unmapped_query_gap_bp_distribution"] == "0;0;0;0;0",
              "Delivered C2 summary has uniform 2-bp overlap and zero query gap", checks)
        read_ends = tsv_rows(delivery / "results" / "c2_read_end_junction_resolution.tsv")
        check(len(read_ends) == 5 and all(row["query_overlap_sequence"] == "AT" for row in read_ends),
              "All five delivered C2 read ends retain AT overlap sequence", checks)
        check(all(row["unmapped_query_gap_bp"] == "0" for row in read_ends),
              "All five delivered C2 read ends retain zero query gap", checks)

        caller = tsv_rows(delivery / "results" / "candidate_manta_delly_concordance.tsv")
        check(len(caller) == 20, "Caller-concordance table has twenty DLD1 comparisons", checks)
        check(sum(int(row["manta_exact_two_end_record_count"]) for row in caller) == 0 and
              sum(int(row["delly_exact_two_end_record_count"]) for row in caller) == 0 and
              sum(int(int(row["manta_exact_two_end_record_count"]) > 0 and int(row["delly_exact_two_end_record_count"]) > 0)
                  for row in caller) == 0,
              "No exact two-end Manta, DELLY, or both-caller match is misreported", checks)
        check(sum(row["caller_concordance_status"] == "anchor_only_or_nonconcordant" for row in caller) == 5,
              "Exactly five anchor-only/nonconcordant records are retained", checks)

        rescue = tsv_rows(delivery / "results" / "read_rescue_candidate_impact.tsv")
        check(len(rescue) == 4 and all(row["rescued_unmapped_mate_qname_overlap"] == "0" for row in rescue),
              "Rescue-impact table has four zero-overlap C1-C4 rows", checks)
        check(all(row["raw_dual_endpoint_qname_count_before_rescue"] == row["raw_dual_endpoint_qname_count_after_rescue"] for row in rescue),
              "Rescue before/after counts agree per candidate", checks)

        clone = tsv_rows(delivery / "results" / "clone_specific_comparison_with_criteria.tsv")
        check(len(clone) == 20, "Clone-comparison table has twenty candidate-by-sample rows", checks)
        c2_s05 = [row for row in clone if row["candidate_label"] == "C2" and row["sample_id"] == "S05_DLD1_c4_dox"]
        check(len(c2_s05) == 1 and c2_s05[0]["raw_dual_endpoint_fragment_count"] == "3" and
              c2_s05[0]["standard_eligible_dual_fragment_count"] == "0",
              "C2 S05 raw and standard-eligible support are correctly separated", checks)

        availability = tsv_rows(delivery / "results" / "data_delivery_availability_manifest.tsv")
        bam = [row for row in availability if row["file_type"] == "BAM"]
        snv = [row for row in availability if row["file_type"] == "project_wide_SNV_indel_VCF"]
        check(len(availability) == 109 and len(bam) == 18 and all(row["file_exists"] == "1" and row["index_exists"] == "1" for row in bam),
              "Manifest retains all 18 existing BAMs with indexes", checks)
        check(sum(int(row["size_bytes"]) for row in bam) == 907440957244,
              "Manifest BAM byte total agrees with report", checks)
        check(len(snv) == 1 and snv[0]["file_exists"] == "0",
              "Manifest explicitly preserves current SNV/indel VCF unavailability", checks)

        for generated in ["customer_question_summary.tsv", "claim_evidence_matrix.tsv"]:
            check((delivery / "results" / generated).is_file(), f"Review-oriented {generated} is supplied", checks)
        check((delivery / "methods" / "source_manifest.json").is_file(), "Source-provenance manifest is supplied", checks)

        forbidden_extensions = {".bam", ".cram", ".fastq", ".fq", ".vcf", ".bcf"}
        disallowed = [str(path.relative_to(delivery)) for path in delivery.rglob("*")
                      if path.is_file() and path.suffix.lower() in forbidden_extensions]
        check(not disallowed, "Supplement does not duplicate large source sequencing/call files", checks)
        check(not (delivery / "P42_Vedita_18_WGS_1023.md").exists(),
              "Supplement does not duplicate the prior full P42 report", checks)
        success, error = True, None
    except Exception as exc:
        success, error = False, str(exc)

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
