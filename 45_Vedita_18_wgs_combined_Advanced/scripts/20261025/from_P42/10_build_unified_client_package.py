#!/usr/bin/env python3
"""Build the one-folder, client-facing P42 delivery after data transfer completes."""
from __future__ import annotations

import argparse
import csv
import json
import shutil
from datetime import date
from pathlib import Path


PROJECT = Path("/home/gao/projects_2026Q4/42_Vedita_18_wgs_advanced_analysis")
CORE = PROJECT / "custom_research_report_20261023"
ASSAY = PROJECT / "custom_research_report_supplementary_20261024_v5"
REQUEST = PROJECT / "custom_research_report_supplementary_20261024_v8"
REPORT = "P42_Vedita_18_WGS_1024.md"


def args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--delivery", type=Path, required=True)
    parser.add_argument("--conversion-audit", type=Path, required=True)
    parser.add_argument("--header-sanitization-audit", type=Path, required=True)
    parser.add_argument("--bam-sanitization-audit", type=Path, required=True)
    return parser.parse_args()


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, str]], fields: list[str]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def copy_file(source: Path, destination: Path) -> None:
    if not source.is_file():
        raise FileNotFoundError(source)
    if destination.exists():
        raise FileExistsError(destination)
    shutil.copy2(source, destination)


def report_text() -> str:
    today = date.today().isoformat()
    return f"""# P42 whole-genome sequencing advanced analysis and P41 Section 7.3 follow-up

**Report Date:** {today}  
**Prepared by:** Zhen Gao, PhD, Principal Bioinformatics Scientist, Athenomics  
**Analysis Platform:** Linux high-performance computing server  
**Reference Genome:** *Homo sapiens* GRCh38  
**Sample context:** DLD1 parental population and clone samples, with and without doxycycline, as named in the follow-up request; tissue of origin was not specified by the client.  
**Delivery status:** Prepared for scientific review and client delivery; not sent to the client.

## 1. Objectives

This single delivery consolidates the complete P42 analysis and the follow-up research requested above P41 report Section 7.3. It provides the full 18-sample short-read whole-genome sequencing (WGS) evidence review, the refined C2 junction analysis, caller concordance, completed unmapped-mate-rescue impact, DLD1 control-versus-doxycycline comparison, an in-silico C2 assay design, and the requested transferable source files.

## 2. Core conclusions

1. P42 evaluated 18 WGS samples across the requested chr11 repeat-region evidence window. The completed inventory retained low-MAPQ, multi-mapping, secondary, supplementary, soft-clipped, and unmapped-mate contexts rather than treating their absence from a filtered candidate list as biological absence.
2. P41 EV0623 is carried forward as P42 C2, not a new event: P41 reported chr11:114,607,192 to chr21:20,648,241, whereas the P42 representation is chr11:114,607,191 to chr21:20,648,242. The one-base difference at each end is a representation difference for the same geometry.
3. C2 is now resolved to five split-read ends from three distinct QNAME labels. All five have a 2-base-pair (bp) `AT` query overlap and a 0-bp query gap at MAPQ 60. This is internally consistent alignment evidence, not independent biological confirmation in a repeat context.
4. Across 20 fixed candidate-by-DLD1-sample comparisons, exact two-end concordance is 0/20 for Manta and 0/20 for DELLY. Five C4 DELLY rows are chr11-only anchors without the chr3 partner and are not called concordant.
5. None of C1-C4 passes the pre-specified high-confidence clone/doxycycline gate: at least three standard-eligible independent QNAME fragments, at least 60x local-depth context for the three-fragment/5% sensitivity target, and sequence-level junction evidence.
6. The C2 assay pair is an in-silico design only. The requested 18 BAM files and indexes, Manta VCFs, and DELLY records are included here. A provenance-qualified project-wide SNV/indel VCF was not identified and is explicitly unavailable.

## 3. Study inputs and sample information

The analysis includes 18 named WGS samples. Mean depth ranges from 8.2x to 19.3x across the complete cohort. The DLD1 follow-up comparison uses S01-S05.

| Sample | Requested role | Treatment label | Use in P41 Section 7.3 follow-up |
| :--- | :---: | :---: | :--- |
| S01_DLD1_pop | Parental population | Not stated as +doxycycline | Comparison reference |
| S02_DLD1_c1_ctrl | Control clone | Control | Comparison reference |
| S03_DLD1_c2_ctrl | Control clone | Control | Comparison reference |
| S04_DLD1_c3_dox | DLD1 clone | +Doxycycline | Case group |
| S05_DLD1_c4_dox | DLD1 clone | +Doxycycline | Case group |

All 18 BAM files and their indexes are included in `genome_wide_data/bam/`. The source BAM total is 907,440,957,244 bytes (907.44 GB; 845.12 GiB). Client BAM headers have local execution paths removed; alignment records and their per-file counts are retained and audited.

## 4. Analysis design and decision criteria

| Research question | Evidence unit and retained context | Decision criterion | Permitted conclusion |
| :--- | :---: | :---: | :--- |
| Does a window contain mapping evidence? | Original alignment records, retaining primary, secondary, supplementary, soft-clipped, multi-mapping, and unmapped-mate context | Evidence counts and filter-impact audit before candidate restriction | Mapping evidence, not an independent molecular count |
| Is C2 read-end structure internally consistent? | Sequence-bearing records grouped by QNAME and read end, retaining original FLAG, CIGAR, MAPQ, and coordinates | Query overlap/gap and overlap sequence across resolved ends | Consistent split-read representation, not a uniquely proven breakpoint |
| Do Manta and DELLY match C1-C4? | All audited records from 18 samples and four caller stages | Both BND ends within 50 bp in either order | Algorithmic concordance from the same WGS data only |
| Is a DLD1 candidate high-confidence? | Fixed-candidate raw and standard-eligible QNAME counts plus local depth | At least 3 standard-eligible QNAMEs, at least 60x context, and sequence-level junction | High-confidence status only if every gate passes |
| Did rescue alter C1-C4? | Recovered unmapped-mate QNAMEs intersected with candidate QNAMEs | Before/after candidate counts and overlap | Measured rescue impact on C1-C4 |

## 5. Methods

P42 evaluated 59 assemblies yielding 90 contigs and 32 alternative placements. From 282 candidate observations, 94 were parseable under the fixed evidence model; four case-only candidates (C1-C4) were retained for detailed review. The raw evidence ledger preserves alignment context and records filter impact rather than allowing a prior caller table to define the data universe.

For the P41 Section 7.3 follow-up, P41 EV0623 was crosswalked to C2. Sequence-bearing alignment records were resolved at `(QNAME, read end)` level using CIGAR-derived query intervals. Manta and DELLY records were compared with each fixed candidate under a two-end, 50-bp tolerance rule. A single nearby anchor, a different partner, or a different SV type remained a retained nonconcordant record. The rescue analysis intersected the completed unmapped-mate ledger with the C1-C4 QNAME sets. The C2 assay is an in-silico primer design and specificity screen; no new assay was performed.

## 6. Results

### 6.1 Complete P42 evidence inventory

The full evidence inventory and assemblies support a repeat-aware interpretation. None of C1-C4 is classified as a high-confidence structural variant. Counts, placements, and candidate-level support are supplied in `results/` and Figures 1-12.

### 6.2 P41 Section 7.3 refinement: C2 read-end junction structure

P41 EV0623 maps to the P42 C2 geometry. Three QNAME labels yielded five resolved split-read ends. Every end has MAPQ 60 for both mapped blocks, a 2-bp `AT` query overlap, and no unaligned query gap. This narrows the technical description beyond the P41 split-read count while preserving the limitation that a repeat-context short-read alignment is not an independent biological confirmation.

### 6.3 Manta and DELLY comparison

| Candidate | Exact Manta two-end matches across S01-S05 | Exact DELLY two-end matches across S01-S05 | Both callers exact | Anchor-only/nonconcordant rows |
| :--- | :---: | :---: | :---: | :---: |
| C1 | 0 | 0 | 0 | 0 |
| C2 | 0 | 0 | 0 | 0 |
| C3 | 0 | 0 | 0 | 0 |
| C4 | 0 | 0 | 0 | 5 |
| Total | 0/20 | 0/20 | 0/20 | 5/20 |

The five C4 nonconcordant rows are recurrent DELLY chr11 insertion anchors without the chr3 partner. They are retained as caller context and are not promoted to evidence for a chr11-chr3 BND.

### 6.4 Unmapped-mate rescue impact

Five unmapped-mate QNAMEs were recovered across the full locus. None overlaps a C1-C4 dual-endpoint QNAME; raw support is unchanged: C1 1 to 1, C2 3 to 3, C3 2 to 2, and C4 1 to 1.

### 6.5 DLD1 control versus doxycycline comparison

| Candidate | Raw support in S04 / S05 | Standard-eligible support in S04 / S05 | Raw support in S01-S03 | Conclusion |
| :--- | :---: | :---: | :---: | :--- |
| C1 | 0 / 1 | 0 / 1 | 0 / 0 / 0 | High-confidence gate not met |
| C2 | 0 / 3 | 0 / 0 | 0 / 0 / 0 | P41 carry-forward hypothesis; high-confidence gate not met |
| C3 | 1 / 1 | 0 / 1 | 0 / 0 / 0 | High-confidence gate not met |
| C4 | 0 / 1 | 0 / 1 | 0 / 0 / 0 | High-confidence gate not met |

For C2, S05 has a raw-support/depth proxy of 3 / 14.760 = 0.203252. It is not a validated allele fraction. C2 standard local depths in S01-S03 are 21.146x, 19.126x, and 23.101x; zero counts at these depths cannot be interpreted as high-sensitivity exclusion at the 60x gate.

### 6.6 In-silico C2 assay design

The selected C2 junction-spanning primer pair is `GTAATCACGGTTACCCAGTGGAT` and `AGCAATAGTTCTCAACAGTGGGT`, with an expected 274-bp product on the virtual junction template. This is an in-silico design and specificity screen only; it does not establish assay performance or confirm C2.

## 7. Figures and supplied results

Figures 1-12 document the complete P42 evidence, assembly, placement, and candidate review. Figure 13 documents the C2 in-silico assay design. Figure 14 summarizes the P41 Section 7.3 client-question follow-up. Every figure is provided in PDF and PNG in `figures/`; `figures/figure_manifest.json` records the 14 formal figures. Machine-readable result tables are in `results/`, including the evidence crosswalk, C2 read-end ledger, caller comparison, rescue impact, clone comparison, assay selection, and client-request coverage.

## 8. Interpretation boundaries and next evidence

### 8.1 What this analysis establishes

The evidence establishes a precise technical characterization of C2, the absence of exact two-end caller concordance under the stated rule, the measured zero C1-C4 rescue increment, and the fact that the fixed candidates do not meet the high-confidence clone/doxycycline gate. It also supplies the specified source-data file categories in a single package.

### 8.2 Claim-to-evidence map

| Claim | Direct evidence | Interpretation boundary |
| :--- | :---: | :--- |
| C2 refines P41 EV0623 | `results/c2_junction_resolution_summary.tsv`; `results/c2_read_end_junction_resolution.tsv` | Same WGS alignment evidence, not independent validation |
| C2 has 2-bp `AT` overlap and 0-bp gap | Five resolved read ends in `results/c2_read_end_junction_resolution.tsv` | Repeat-context placement remains ambiguous |
| No exact Manta or DELLY support | `results/candidate_manta_delly_concordance.tsv` | Caller concordance would still be algorithmic evidence |
| Rescue did not increase C1-C4 support | `results/read_rescue_candidate_impact.tsv` | Does not exclude rescue value elsewhere |
| Clone/doxycycline gate is not met | `results/clone_specific_comparison_with_criteria.tsv` | Low local depth limits negative sensitivity |
| C2 assay is design-only | `results/assay_design_selection.tsv`; `methods/virtual_junction_templates.fasta` | No independent assay result is included |

### 8.3 Limitations

This is a deeper analysis of the same short-read WGS data, not an independent biological replicate. A split-read overlap supports a mapping representation but cannot alone prove a unique breakpoint in a repetitive locus. No independent breakpoint-spanning assay, new sequencing, raw fluorescence in situ hybridization input, or scored cytogenetic result was supplied for this package. The unavailable project-wide SNV/indel VCF is not substituted by Manta candidate small-indel output.

### 8.4 Next evidence

The appropriate next evidentiary step for C2 is experimental assessment of the supplied junction-spanning assay design and, if available, integration of independently generated cytogenetic scoring with pre-specified interpretation rules. Any later analysis should retain the raw-alignment and filter-impact audit used here.

## 9. Delivery contents

| Directory or file | Contents | Status |
| :--- | :---: | :--- |
| `P42_Vedita_18_WGS_1024.md` | This integrated client-facing report | Included |
| `figures/` | 14 formal figures, each as PDF and PNG | Included |
| `results/` | Machine-readable analysis and follow-up tables | Included |
| `methods/` | Client-safe methods manifest and virtual junction template | Included |
| `genome_wide_data/bam/` | 18 BAM files and 18 BAM indexes; alignment records preserved and local execution paths removed from headers | Included |
| `genome_wide_data/manta/` | Filtered SV, candidate SV, and candidate small-indel VCFs; variant records preserved and local execution paths removed from headers | Included |
| `genome_wide_data/delly/` | Discovery and region-genotyped VCFs, converted from source BCF with indexes; variant records preserved and local execution paths removed from headers | Included |
| `results/data_inventory.tsv` | Client-relative inventory of all supplied data | Included |
| `methods/vcf_header_sanitization_audit.tsv` | Per-file variant-body hash comparison and removed-header-path count | Included |
| `methods/bam_header_sanitization_audit.tsv` | Per-file alignment-count comparison and removed-header-path count | Included |
| Project-wide SNV/indel VCF | No provenance-qualified source identified | Unavailable; not represented as present |
"""


def main() -> None:
    ns = args()
    delivery = ns.delivery.resolve()
    if delivery.name != "custom_research_report_20261024":
        raise ValueError("Delivery directory must be custom_research_report_20261024")
    required_dirs = [delivery / name for name in ("figures", "results", "methods", "genome_wide_data")]
    if not all(path.is_dir() for path in required_dirs):
        raise FileNotFoundError("Expected skeleton and copied data are incomplete")
    for name in (REPORT, "figures/figure_manifest.json", "results/data_inventory.tsv"):
        if (delivery / name).exists():
            raise FileExistsError(f"Refusing to overwrite {delivery / name}")

    figure_sources = sorted((CORE / "figures").glob("P42_Fig*.*"))
    figure_sources += sorted((ASSAY / "figures").glob("P42_Fig13*.*"))
    figure_sources += sorted((REQUEST / "figures").glob("P42_Fig14*.*"))
    figures = delivery / "figures"
    for source in figure_sources:
        copy_file(source, figures / source.name)
    pngs = sorted(figures.glob("*.png"))
    pdfs = sorted(figures.glob("*.pdf"))
    if len(pngs) != 14 or len(pdfs) != 14:
        raise ValueError(f"Expected 14 PNG and 14 PDF figures, got {len(pngs)} and {len(pdfs)}")
    (figures / "figure_manifest.json").write_text(json.dumps({"figure_count": 14, "png": [p.name for p in pngs], "pdf": [p.name for p in pdfs]}, indent=2) + "\n")

    tables = [
        (CORE / "results" / "contig_placement_summary.tsv", "contig_placement_summary.tsv"),
        (CORE / "results" / "fixed_candidate_dld1_matrix.tsv", "fixed_candidate_dld1_matrix.tsv"),
        (CORE / "results" / "fixed_candidate_summary.tsv", "fixed_candidate_summary.tsv"),
        (CORE / "results" / "raw_evidence_match_audit.tsv.gz", "raw_evidence_match_audit.tsv.gz"),
        (ASSAY / "results" / "assay_design_selection.tsv", "assay_design_selection.tsv"),
        (ASSAY / "results" / "candidate_template_audit.tsv", "candidate_template_audit.tsv"),
        (ASSAY / "results" / "external_input_status.tsv", "external_input_status.tsv"),
        (ASSAY / "results" / "primer_specificity.tsv", "primer_specificity.tsv"),
        (ASSAY / "results" / "raw_support_alignment_audit.tsv", "raw_support_alignment_audit.tsv"),
        (REQUEST / "results" / "c2_junction_resolution_summary.tsv", "c2_junction_resolution_summary.tsv"),
        (REQUEST / "results" / "c2_read_end_junction_resolution.tsv", "c2_read_end_junction_resolution.tsv"),
        (REQUEST / "results" / "candidate_manta_delly_concordance.tsv", "candidate_manta_delly_concordance.tsv"),
        (REQUEST / "results" / "claim_evidence_matrix.tsv", "claim_evidence_matrix.tsv"),
        (REQUEST / "results" / "clone_specific_comparison_with_criteria.tsv", "clone_specific_comparison_with_criteria.tsv"),
        (REQUEST / "results" / "customer_question_summary.tsv", "customer_question_summary.tsv"),
        (REQUEST / "results" / "customer_request_coverage.tsv", "customer_request_coverage.tsv"),
        (REQUEST / "results" / "read_rescue_candidate_impact.tsv", "read_rescue_candidate_impact.tsv"),
    ]
    for source, name in tables:
        copy_file(source, delivery / "results" / name)
    copy_file(ASSAY / "methods" / "virtual_junction_templates.fasta", delivery / "methods" / "virtual_junction_templates.fasta")
    copy_file(ns.conversion_audit, delivery / "methods" / "delly_conversion_audit.tsv")
    copy_file(ns.header_sanitization_audit, delivery / "methods" / "vcf_header_sanitization_audit.tsv")
    copy_file(ns.bam_sanitization_audit, delivery / "methods" / "bam_header_sanitization_audit.tsv")

    inventory_rows: list[dict[str, str]] = []
    for category, directory, pattern, index_suffix in [
        ("BAM", "bam", "*.bam", ".bai"),
        ("Manta_filtered_SV", "manta/filtered_sv", "*.vcf.gz", ".tbi"),
        ("Manta_candidate_SV", "manta/candidate_sv", "*.vcf.gz", ".tbi"),
        ("Manta_candidate_small_indel", "manta/candidate_small_indel", "*.vcf.gz", ".tbi"),
        ("DELLY_discovery_VCF", "delly/discovery", "*.vcf.gz", ".tbi"),
        ("DELLY_region_genotyped_VCF", "delly/region_genotyped", "*.vcf.gz", ".tbi"),
    ]:
        for path in sorted((delivery / "genome_wide_data" / directory).glob(pattern)):
            inventory_rows.append({"category": category, "client_relative_path": str(path.relative_to(delivery)), "bytes": str(path.stat().st_size), "index_client_relative_path": str(Path(str(path) + index_suffix).relative_to(delivery)), "availability": "included"})
    inventory_rows.append({"category": "project_wide_SNV_indel_VCF", "client_relative_path": "", "bytes": "0", "index_client_relative_path": "", "availability": "unavailable: no provenance-qualified source identified"})
    write_tsv(delivery / "results" / "data_inventory.tsv", inventory_rows, ["category", "client_relative_path", "bytes", "index_client_relative_path", "availability"])
    manifest = {"package": delivery.name, "report": REPORT, "figure_count": 14, "data_categories": sorted({r["category"] for r in inventory_rows}), "source_paths_redacted": True}
    (delivery / "methods" / "p42_methods_and_decision_criteria.json").write_text(json.dumps(manifest, indent=2) + "\n")
    (delivery / REPORT).write_text(report_text())


if __name__ == "__main__":
    main()
