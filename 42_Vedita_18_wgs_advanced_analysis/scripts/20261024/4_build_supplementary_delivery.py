#!/usr/bin/env python3
"""Build the 2026-10-24 incremental delivery for the P41 Section 7.3 follow-up.

This script deliberately copies only new, date-stamped analysis outputs.  It never
rewrites the 2026-10-23 delivery and refuses to run if the supplementary directory
already exists.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import shutil
from datetime import date
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project", type=Path, required=True)
    parser.add_argument("--out-dir", type=Path, required=True)
    parser.add_argument("--delivery", type=Path, required=True)
    return parser.parse_args()


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, str]], fields: list[str]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def selection_summary(rows: list[dict[str, str]]) -> dict[str, dict[str, str]]:
    return {row["candidate_label"]: row for row in rows}


def main() -> None:
    args = parse_args()
    out_dir = args.out_dir.resolve()
    delivery = args.delivery.resolve()
    if delivery.exists():
        raise FileExistsError(f"Refusing to overwrite existing delivery directory: {delivery}")

    required = [
        out_dir / "assay_design_selection.tsv",
        out_dir / "primer_specificity.tsv",
        out_dir / "candidate_template_audit.tsv",
        out_dir / "raw_support_alignment_audit.tsv",
        out_dir / "virtual_junction_templates.fasta",
        out_dir / "figures" / "P42_Fig13_P41_7_3_assay_followup.png",
        out_dir / "figures" / "P42_Fig13_P41_7_3_assay_followup.pdf",
    ]
    missing = [str(path) for path in required if not path.is_file()]
    if missing:
        raise FileNotFoundError("Required follow-up output(s) missing:\n" + "\n".join(missing))

    delivery.mkdir()
    results = delivery / "results"
    figures = delivery / "figures"
    methods = delivery / "methods"
    results.mkdir()
    figures.mkdir()
    methods.mkdir()

    copies = {
        out_dir / "assay_design_selection.tsv": results / "assay_design_selection.tsv",
        out_dir / "primer_specificity.tsv": results / "primer_specificity.tsv",
        out_dir / "candidate_template_audit.tsv": results / "candidate_template_audit.tsv",
        out_dir / "raw_support_alignment_audit.tsv": results / "raw_support_alignment_audit.tsv",
        out_dir / "virtual_junction_templates.fasta": methods / "virtual_junction_templates.fasta",
        out_dir / "figures" / "P42_Fig13_P41_7_3_assay_followup.png": figures / "P42_Fig13_P41_7_3_assay_followup.png",
        out_dir / "figures" / "P42_Fig13_P41_7_3_assay_followup.pdf": figures / "P42_Fig13_P41_7_3_assay_followup.pdf",
    }
    for source, destination in copies.items():
        shutil.copy2(source, destination)

    status_rows = [
        {
            "required_external_input": "Independent breakpoint-spanning PCR result from DLD1 clone 4 ± doxycycline genomic DNA",
            "status": "not_provided",
            "reason": "No new wet-lab result was available to this computational follow-up.",
            "effect_on_interpretation": "EV0623/C2 remains a candidate mapping event, not a confirmed clone-specific rearrangement.",
        },
        {
            "required_external_input": "Local WGS depth at or above 60× for sensitivity testing of low-fraction C2 support",
            "status": "not_provided",
            "reason": "No new sequencing data were available.",
            "effect_on_interpretation": "Current P42 short-read evidence cannot exclude a low-fraction C2 event at the specified sensitivity target.",
        },
        {
            "required_external_input": "Raw FISH images and per-cell scoring table",
            "status": "not_provided",
            "reason": "No raw image data or scoring table were available.",
            "effect_on_interpretation": "No association between the candidate locus and the observed cytogenetic subpopulation was tested.",
        },
    ]
    write_tsv(
        results / "external_input_status.tsv",
        status_rows,
        ["required_external_input", "status", "reason", "effect_on_interpretation"],
    )

    selected = selection_summary(read_tsv(out_dir / "assay_design_selection.tsv"))
    c2 = selected.get("C2", {})
    c2_status = c2.get("in_silico_design_status", "not_available")
    c2_pair = c2.get("selected_unchecked_pair_rank", "none") or "none"
    c2_note = c2.get("interpretation", "No C2 assay selection record was produced.")
    c2_raw_support = c2.get("raw_dual_endpoint_qname_count", "0")
    c2_concordant = c2.get("raw_junction_concordant_qname_count", "0")
    c2_endpoints = c2.get("fixed_endpoints", "chr11:114,607,191 / chr21:20,648,242")
    c2_left = c2.get("left_primer", "not selected") or "not selected"
    c2_right = c2.get("right_primer", "not selected") or "not selected"
    c2_product = c2.get("expected_junction_product_bp", "not selected") or "not selected"
    report = f"""# P42 supplementary follow-up: P41 Section 7.3 candidate-event assay design

**Client/project:** Vedita-18 WGS advanced analysis (P42)  
**Date:** {date.today().isoformat()}  
**Species:** *Homo sapiens* (GRCh38)  
**Tissue/cell context:** DLD1 clone 4 ± doxycycline and P41 repeat-region follow-up; no additional biological material was analysed in this computational supplement.  
**Prepared by:** Zhen Gao, PhD, Principal Bioinformatics Scientist, Athenomics

## 1. Objective

This supplementary analysis is an incremental follow-up to Section 7.3 of the P41 Chr11 repeat-region WGS report. P41 identified EV0623, a low-support chr11–chr21 split-read candidate in DLD1 clone 4 + doxycycline. P42 reconciled that event with candidate C2 (one-base endpoint representation difference at each breakpoint end) but did not establish a high-confidence clone- or doxycycline-specific rearrangement. The present work converts the C2 geometry into an orientation-aware, reference-derived in-silico assay-design record and performs an in-silico specificity screen. It does not repeat the completed cohort-wide P42 analysis.

## 2. Inputs and data scope

The input evidence was restricted to the frozen P42 candidate ledger, normalized SV-event table, and the existing S05 aligned-read record used to audit raw support. The same GRCh38 reference assembly was used for virtual junction construction and primer-specificity assessment. No new FASTQ, BAM, PCR product, sequencing run, or raw FISH image was supplied.

## 3. Key findings

### What P42 adds beyond P41 Section 7.3

1. **From a candidate observation to an exact crosswalk.** P41 reported EV0623 as a DLD1 clone 4 + doxycycline candidate supported by 3 split reads and 0 discordant pairs. P42 reconciled that observation with C2 at **{c2_endpoints}**. The one-base representation differences relative to P41 are breakpoint representation differences, not a second event.
2. **From caller support to retained raw-sequence evidence.** Instead of treating supplementary alignments as independent molecules, P42 preserved their QNAME, FLAG, CIGAR, MAPQ, and SA attributes and tested the recovered read sequence against alternative virtual-junction orientations. The final raw-sequence audit finds **{c2_concordant}/{c2_raw_support}** C2 dual-endpoint QNAMEs with exact junction-spanning 24-mer concordance under the selected template orientation.
3. **From a locus hypothesis to a coordinate- and strand-constrained assay.** P42 generated candidate junction PCR pairs, rejected identical-primer and excessive heterodimer combinations, and required each selected primer to have one exact full-length GRCh38 hit at its expected chromosome, coordinate interval, and strand. C2 selection status: **{c2_status}**; selected pair rank: **{c2_pair}**; expected product: **{c2_product} bp**; primers: **{c2_left}** / **{c2_right}**. Selection note: {c2_note}
4. **Stronger technical resolution, unchanged biological caution.** These analyses make C2 more precisely defined and directly testable than the P41 observation, but do not establish that an amplicon exists in independent clone 4 + doxycycline DNA, that C2 is clone-specific, or that C2 explains the P41 copy-number or cytogenetic observations.

## 4. Results overview

Figure 1 summarizes the evidence chain from the P41 EV0623 observation through the P42 C2 crosswalk, raw-alignment/template audit, in-silico assay screen, and the external observations still required for biological validation. The machine-readable candidate, raw-alignment, primer-specificity, and selection tables are supplied in `results/`.

![Figure 1: P41 Section 7.3 assay follow-up]({figures.name}/P42_Fig13_P41_7_3_assay_followup.png)

## 5. Candidate-level result

The C2 virtual junction was constructed from the normalized reciprocal breakend coordinates, a small coordinate-representation grid, and both possible remote reference orientations. The selected orientation was the one with maximum exact support in recovered raw split-read sequences. The audit preserves original alignment attributes (QNAME, FLAG, coordinate, CIGAR, MAPQ, and SA tag where present) rather than treating supplementary alignments as independent molecules. C1, C3, and C4 were retained as contextual P42 candidates and were evaluated by the same design pipeline, but are not reclassified by this supplement.

## 6. Methods

For each candidate, virtual breakpoint templates were assembled from GRCh38 flanking sequences across a coordinate grid and both remote reference orientations. The C2 raw-support alignments were retrieved from the source BAM in local breakpoint windows and compared to junction-spanning template k-mers; the orientation with maximum raw-sequence concordance was retained. Primer3-py generated short amplicon candidates around each selected virtual junction. Candidate pairs with identical primer sequences or excessive predicted heterodimerization were rejected. Each remaining primer was screened against a locally built GRCh38 BLAST database using `blastn-short`; only exact full-length single-locus hits at the expected chromosome, coordinate interval, and strand were eligible for the in-silico selection status. This is a computational specificity screen, not laboratory validation.

## 7. Interpretation

The P41-to-P42 event reconciliation is technically coherent, and the supplied assay tables make the C2 hypothesis testable with independent material. However, low raw read support, repeat-context mapping ambiguity, and the absence of a high-depth sensitivity dataset mean that the result remains hypothesis-generating. No conclusion is made about causality for the P41 copy-number step or FISH-defined subpopulation.

## 8. Limitations and next work

No PCR, new sequencing, or raw FISH quantification was performed in this computational supplementary analysis. No independent breakpoint-spanning PCR result, new ≥60× local WGS data, or raw FISH image/scoring data was available. The required external inputs and their impact are listed in `results/external_input_status.tsv`.

The next evidence-bearing steps are: (1) test the selected C2 primer pair, if any, against independent DLD1 clone 4 ± doxycycline genomic DNA with suitable controls; (2) if sensitivity remains limiting, assess new local short-read WGS depth against the pre-specified ≥60× target; and (3) evaluate raw FISH images and per-cell scores before inferring a relationship to the cytogenetic observation. These steps require new experimental material or source data and were not performed here.

## 9. Deliverables

This incremental package contains one publication-style figure (PNG and PDF), candidate/template and raw-alignment audit tables, primer-specificity and assay-selection tables, a virtual-junction FASTA file, this report, and a machine-readable list of external inputs still required. It intentionally contains no duplicate P42 cohort-wide results and does not alter the P41 or 2026-10-23 P42 delivery directories.

---

Prepared for scientific review; not sent to the client.
"""
    report_path = delivery / "P42_P41_7_3_Assay_Followup_1024.md"
    report_path.write_text(report)

    manifest = {
        "analysis_date": date.today().isoformat(),
        "scope": "Incremental computational follow-up to P41 report Section 7.3; does not repeat or replace P42 2026-10-23 delivery.",
        "source_out_dir": str(out_dir),
        "delivery": str(delivery),
        "source_sha256": {str(source.relative_to(out_dir)): sha256(source) for source in copies},
        "candidate_count": len(selected),
        "c2_selection_status": c2_status,
        "c2_selected_unchecked_pair_rank": c2_pair,
    }
    (methods / "assay_design_manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")


if __name__ == "__main__":
    main()
