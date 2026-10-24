#!/usr/bin/env python3
"""Build an append-only P42 client-question supplementary package."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import shutil
from datetime import date
from pathlib import Path


REPORT_NAME = "P42_Client_Request_Completion_1024.md"
FIGURE_STEM = "P42_Fig14_Client_Request_Completion"


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
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


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def int_or_zero(value: str) -> int:
    return int(value) if value not in {"", "."} else 0


def source_manifest(rows: list[dict[str, str]]) -> dict[str, dict[str, int]]:
    result: dict[str, dict[str, int]] = {}
    for file_type in sorted({row["file_type"] for row in rows}):
        subset = [row for row in rows if row["file_type"] == file_type]
        result[file_type] = {
            "entries": len(subset),
            "existing_files": sum(int_or_zero(row["file_exists"]) for row in subset),
            "indexes_present": sum(int_or_zero(row["index_exists"]) for row in subset),
            "size_bytes": sum(int_or_zero(row["size_bytes"]) for row in subset),
        }
    return result


def make_report(manifest: dict[str, dict[str, int]]) -> str:
    bam = manifest["BAM"]
    manta_filtered = manifest["Manta_filtered_SV"]
    manta_candidate_sv = manifest["Manta_candidate_SV"]
    manta_candidate_indel = manifest["Manta_candidate_small_indel"]
    delly_discovery = manifest["DELLY_discovery_BCF"]
    delly_genotyped = manifest["DELLY_region_genotyped_BCF"]
    bam_gb = bam["size_bytes"] / 1e9
    bam_gib = bam["size_bytes"] / (2 ** 30)
    return f"""# P42 supplementary completion: client questions on P41 Section 7.3

**Report Date:** {date.today().isoformat()}  
**Prepared by:** Zhen Gao, PhD, Principal Bioinformatics Scientist, Athenomics  
**Analysis Platform:** Linux HPC server  
**Species:** *Homo sapiens* (GRCh38)  
**Tissue / Cell:** DLD1 parental population and DLD1 clone samples, with and without doxycycline, as named in the client request  

## 1. Objectives

This incremental supplement answers the follow-up questions directed to Section 7.3 of the P41 chr11 repeat-region whole-genome sequencing (WGS) report. It uses existing GRCh38-aligned records only; no re-alignment or new sequencing was performed. The new work: (i) resolves the available P42 C2 split-read alignment representation at read-end level, (ii) records whether Manta and DELLY provide exact two-end breakend (BND) concordance, (iii) quantifies the effect of the completed unmapped-mate rescue on the C1-C4 candidates, (iv) applies the pre-specified DLD1 parental/control versus +doxycycline comparison criteria, and (v) inventories the requested genome-wide source files without copying the large source data into this supplementary package.

## 2. Core Conclusions

1. **C2 is now more precisely characterized than in P41.** The P41 EV0623 geometry crosswalks to P42 C2 at chr11:114,607,191 and chr21:20,648,242. Across three distinct SAM read-name (QNAME) labels and five resolved split-read ends, every read-end representation has a 2-base-pair (bp) `AT` query overlap and a 0-bp unaligned query gap. This is stronger technical detail than the P41 split-read count, but it remains repeat-context alignment evidence rather than independent biological confirmation.
2. **None of C1-C4 has exact two-end Manta or DELLY concordance in the audited DLD1 records.** Across 20 candidate-by-sample comparisons, exact two-end matches were 0/20 for Manta, 0/20 for DELLY, and 0/20 for both callers. Five C4 rows have a recurrent DELLY chr11 insertion anchor without the chr3 partner; they are deliberately classified as anchor-only/nonconcordant rather than as support for C4.
3. **No C1-C4 event reaches the P42 high-confidence clone- or doxycycline-specific gate.** C2 has three raw dual-endpoint QNAMEs only in S05, but zero standard-eligible fragments and local depth below the pre-specified 60x sensitivity context. The observed raw-support/depth ratio is a descriptive proxy, not a validated allele fraction.
4. **Read rescue did not increase C1-C4 support.** Five unmapped-mate QNAMEs were retained across the locus, yet none overlap a C1-C4 dual-endpoint QNAME; each candidate therefore has the same raw count before and after rescue.
5. **The requested source-data inventory is ready for administrative transfer planning, with one explicit gap.** All 18 BAM files and all 18 BAM indexes are present (907.44 GB; 845.12 GiB). Existing Manta and DELLY source files are inventoried in this package. No provenance-qualified project-wide single-nucleotide variant (SNV)/small indel Variant Call Format (VCF) source was identified, so no such file is represented as available.

Section 3 organizes the measured findings by theme. Section 8.2 gives the evidence path and interpretation boundary for each conclusion; Section 8.3 states what the data do not establish, and Section 8.4 separates scientific next evidence from administrative actions.

## 3. Key Findings

1. **P41-to-P42 precision:** P42 retains the P41 EV0623 observation and resolves its correspondence to one fixed C2 geometry rather than treating the one-base representation difference at each end as a second event.
2. **Junction structure:** C2 read ends carry a consistent `AT` overlap in the split-alignment representation, with no observed unaligned query bases between the two mapped blocks.
3. **Caller comparison:** Manta/DELLY were assessed under a two-end, 50-bp tolerance rule, so a nearby single-anchor insertion is not silently promoted to a trans-chromosomal BND match.
4. **Clone comparison:** The client-requested five DLD1 samples were evaluated against the same candidate set, support rule, depth context, and sequence-level requirement; low local depth prevents a negative call from being interpreted as exclusion at the specified sensitivity target.
5. **Delivery readiness:** The manifest distinguishes present source files from source files that have not been demonstrated, and it records size and index status instead of duplicating approximately 0.91 TB of BAM data.

## 4. Sample Information

| Sample | DLD1 role in the requested comparison | Treatment label in client request | Role in this supplement |
| :--- | :---: | :---: | :---: |
| S01_DLD1_pop | Parental population | not stated as +doxycycline | Comparison reference |
| S02_DLD1_c1_ctrl | Control clone | control | Comparison reference |
| S03_DLD1_c2_ctrl | Control clone | control | Comparison reference |
| S04_DLD1_c3_dox | DLD1 clone | +doxycycline | Case group |
| S05_DLD1_c4_dox | DLD1 clone | +doxycycline | Case group |

## 5. Analysis Rationale and Decision Criteria

| Question | Observation / evidence unit | Retention and counting rule | Decision criterion | Permitted interpretation |
| :--- | :---: | :---: | :---: | :--- |
| Is C2 internally consistent at read-end level? | Sequence-bearing primary and supplementary alignment records grouped by `(QNAME, read end)` | Preserve original FLAG, CIGAR, mapping quality (MAPQ), coordinates, and query intervals; QNAMEs count observations but are not biological replicates | chr11/chr21 query overlap or gap and overlap sequence across all resolved ends | A split-read alignment representation consistent with microhomology; not a uniquely proven breakpoint |
| Do Manta and DELLY concord on C1-C4? | All audited records from 72 source inputs (18 samples x 4 caller stages) | Retain every source record; do not convert a missing record into evidence of biological absence | Both BND ends must match a fixed candidate within 50 bp in either order; a single nearby endpoint or different partner/type is anchor-only/nonconcordant | Algorithmic concordance from the same WGS dataset, not independent validation |
| Is a DLD1 +doxycycline candidate high-confidence? | Fixed-candidate dual-endpoint fragment counts and local depth in S01-S05 | Raw and standard-eligible counts remain separate; repeat-context support/depth is not treated as an allele-fraction estimator | At least 3 standard-eligible independent QNAME fragments, local depth at least 60x for the 3-fragment/5% sensitivity context, and a sequence-level junction | High-confidence status only if all gates pass; no direct allele-fraction cutoff is applied |
| Did unmapped-mate rescue change C1-C4? | Recovered unmapped-mate QNAMEs intersected with each candidate’s raw dual-endpoint QNAMEs | Retain the five recovered mates and report each overlap, including zero | Before/after raw QNAME count and intersection size for each candidate | Measured rescue impact on C1-C4; no overlap does not mean rescue has no value elsewhere |
| Which requested data types are available? | Source file, index, existence flag, and byte size | Inventory is read-only; large BAM/VCF/BCF files are not copied into this scientific supplement | Availability requires a present, provenance-qualified source path; candidate-stage Manta files are not renamed as “unfiltered” | A transfer-planning manifest, with unavailable categories stated explicitly |

## 6. Methods

The work began with the completed P42 whole-window evidence inventory and local assembly/read-rescue outputs, retaining mapping-quality-zero, multi-mapping, secondary, supplementary, soft-clipped, and unmapped-mate contexts where present. This supplement did not re-align the existing BAM files.

For C2, all available sequence-bearing alignment records from the prior raw-support audit were read at the level of individual read ends. Query coordinates were derived from CIGAR strings, then the chr11 and chr21 blocks were compared for overlap and gap. This retains primary and supplementary records as structural context while avoiding a claim that multiple records are independent molecules.

For caller concordance, 72 audited Manta and DELLY source files were compared with each C1-C4 geometry in each DLD1 sample. An exact match required both BND endpoints to fall within 50 bp in either partner order. A nearby single-end anchor or a record without the fixed candidate’s partner was retained and labelled nonconcordant. For the clone comparison, raw dual-endpoint counts, standard-eligible counts, and standard local depth were read from the completed P42 fixed-candidate matrix. The delivery manifest records only source availability and size; it does not alter, copy, or send any source dataset.

## 7. Results

### 7.1 P41 Section 7.3 to P42 C2 refinement

P41 EV0623 was reported at chr11:114,607,192 to chr21:20,648,241 with three split reads and zero discordant pairs in S05. P42 C2 is represented as chr11:114,607,191 to chr21:20,648,242. The one-base difference at each endpoint is a breakpoint-representation difference in the crosswalk of the same candidate geometry. P42 previously established an in-silico C2 assay design record; the new result here is the independent read-end accounting of the existing alignment evidence.

### 7.2 C2 read-end junction detail

| C2 measurement | Observed value | Interpretation boundary |
| :--- | :---: | :--- |
| Distinct QNAMEs | 3 | Counted as retained observations, not independent biological validation |
| Resolved split-read ends | 5 | All are from S05 and retain their alignment context |
| chr11 and chr21 MAPQ | 60 for each resolved block | High alignment confidence does not resolve repeat-context placement ambiguity by itself |
| Query overlap | 2 bp in 5/5 read ends | Overlap sequence is `AT` in every resolved end |
| Unaligned query gap | 0 bp in 5/5 read ends | No intervening unaligned query sequence was observed |
| Allowed conclusion | Alignment representation consistent with 2-bp microhomology | Not an independently verified, unique biological breakpoint |

### 7.3 Manta and DELLY concordance

| Candidate | Exact Manta two-end matches across S01-S05 | Exact DELLY two-end matches across S01-S05 | Both callers exact | Anchor-only / nonconcordant rows | Interpretation |
| :--- | :---: | :---: | :---: | :---: | :--- |
| C1 | 0 | 0 | 0 | 0 | No exact caller match in the audited sources |
| C2 | 0 | 0 | 0 | 0 | No exact caller match in the audited sources |
| C3 | 0 | 0 | 0 | 0 | No exact caller match in the audited sources |
| C4 | 0 | 0 | 0 | 5 | Recurrent DELLY chr11 insertion anchor lacks the chr3 partner and is not C4 concordance |
| Total | 0/20 | 0/20 | 0/20 | 5/20 | No C1-C4 geometry is called exactly by Manta, DELLY, or both under the two-end rule |

### 7.4 Read-rescue impact

| Candidate | Raw dual-endpoint QNAMEs before rescue | Rescued-mate QNAME overlap | Raw dual-endpoint QNAMEs after rescue | Result |
| :--- | :---: | :---: | :---: | :--- |
| C1 | 1 | 0 | 1 | Unchanged |
| C2 | 3 | 0 | 3 | Unchanged |
| C3 | 2 | 0 | 2 | Unchanged |
| C4 | 1 | 0 | 1 | Unchanged |

Five unmapped-mate QNAMEs were recovered across the full locus. None intersected the four candidate-specific dual-endpoint QNAME sets above.

### 7.5 DLD1 parental/control versus +doxycycline comparison

| Candidate | Raw support in S04 / S05 | Standard-eligible support in S04 / S05 | Raw support in S01-S03 | Applicable conclusion |
| :--- | :---: | :---: | :---: | :--- |
| C1 | 0 / 1 | 0 / 1 | 0 / 0 / 0 | Low-support S05 candidate; high-confidence gate not met |
| C2 | 0 / 3 | 0 / 0 | 0 / 0 / 0 | P41 carry-forward candidate in S05; high-confidence gate not met |
| C3 | 1 / 1 | 0 / 1 | 0 / 0 / 0 | Low-support candidates in both +doxycycline samples; high-confidence gate not met |
| C4 | 0 / 1 | 0 / 1 | 0 / 0 / 0 | Low-support S05 candidate; high-confidence gate not met |

For C2, S05 has raw support/depth proxy 3 / 14.760 = 0.203252; S04 has 0 / 9.627. This ratio is displayed only as a raw-support/depth proxy because repeat-context evidence and standard depth do not form a validated allele-fraction denominator. The three comparison references have C2 standard local depth 21.146x, 19.126x, and 23.101x, respectively, below the 60x pre-specified sensitivity context; their zero raw counts must not be read as high-sensitivity exclusion.

### 7.6 Requested genome-wide source-file availability

| Requested data category | Source entries present | Index status recorded | Aggregate source size | Availability statement |
| :--- | :---: | :---: | :---: | :--- |
| BAM | {bam['existing_files']}/{bam['entries']} | {bam['indexes_present']}/{bam['entries']} indexes present | {bam['size_bytes']:,} bytes ({bam_gb:.2f} GB; {bam_gib:.2f} GiB) | Available as source data; not duplicated here |
| Manta filtered SV VCF | {manta_filtered['existing_files']}/{manta_filtered['entries']} | {manta_filtered['indexes_present']}/{manta_filtered['entries']} indexes present | {manta_filtered['size_bytes']:,} bytes | Available filtered source |
| Manta candidate SV VCF | {manta_candidate_sv['existing_files']}/{manta_candidate_sv['entries']} | {manta_candidate_sv['indexes_present']}/{manta_candidate_sv['entries']} index entries recorded | {manta_candidate_sv['size_bytes']:,} bytes | Candidate-stage source; not relabelled as “unfiltered” |
| Manta candidate small-indel VCF | {manta_candidate_indel['existing_files']}/{manta_candidate_indel['entries']} | {manta_candidate_indel['indexes_present']}/{manta_candidate_indel['entries']} indexes present | {manta_candidate_indel['size_bytes']:,} bytes | Candidate small-indel source; not a project-wide SNV/indel callset |
| DELLY discovery Binary Call Format (BCF) | {delly_discovery['existing_files']}/{delly_discovery['entries']} | {delly_discovery['indexes_present']}/{delly_discovery['entries']} index entries recorded | {delly_discovery['size_bytes']:,} bytes | Discovery-stage source |
| DELLY region-genotyped BCF | {delly_genotyped['existing_files']}/{delly_genotyped['entries']} | {delly_genotyped['indexes_present']}/{delly_genotyped['entries']} index entries recorded | {delly_genotyped['size_bytes']:,} bytes | Region-genotyped source |
| Project-wide SNV/indel VCF | 0/1 | 0/1 | 0 bytes | No provenance-qualified source was identified for the current package |

## 8. Conclusions, Limitations and Next Steps

### 8.1 Candidate-level conclusion table

| Candidate | Current status | Evidence supporting the status | What is not established |
| :--- | :---: | :---: | :--- |
| C1 | Low-support S05 candidate | One raw and one standard-eligible dual-endpoint QNAME in S05 | Clone/doxycycline specificity or a confirmed junction |
| C2 | Highest-priority P41 EV0623 carry-forward hypothesis | Three S05 raw QNAMEs, five resolved read ends, uniform 2-bp `AT` overlap, and an in-silico assay design record | Independent breakpoint confirmation, clone specificity, or causal link to cytogenetic observations |
| C3 | Low-support candidate in both +doxycycline samples | One raw QNAME in each +doxycycline sample; one standard-eligible QNAME in S05 | High-confidence event status or a sequence-level junction |
| C4 | Low-support S05 candidate with nonconcordant DELLY anchor context | One raw and one standard-eligible QNAME in S05; recurrent DELLY chr11 anchor-only records | chr11-chr3 BND concordance or clone/doxycycline specificity |

### 8.2 Practical reading with evidence paths

| Core conclusion | Why this is the appropriate reading | Evidence location |
| :--- | :---: | :--- |
| C2 gained base-pair detail but not independent confirmation | Five alignment read ends from three QNAMEs all show `AT` 2-bp overlap and 0-bp gap; repeat-context mapping and the single WGS data source define the boundary | Section 7.2; `results/c2_read_end_junction_resolution.tsv`; `results/c2_junction_resolution_summary.tsv` |
| Exact Manta/DELLY caller concordance is absent | The two-end 50-bp rule prevents a recurrent single chr11 insertion anchor from being conflated with a chr11-chr3 BND | Section 7.3; `results/candidate_manta_delly_concordance.tsv` |
| C1-C4 remain below the clone-specific gate | Raw support, standard-eligible support, local depth, and sequence-level criteria are recorded separately; C2 S05 has zero standard-eligible support | Section 5 and Section 7.5; `results/clone_specific_comparison_with_criteria.tsv` |
| Rescue produced a measured zero increment for C1-C4 | The five recovered unmapped mates are retained in the parent rescue ledger and their candidate intersection is explicitly counted as zero | Section 7.4; `results/read_rescue_candidate_impact.tsv` |
| Transfer planning can start from the manifest, with one gap | The inventory names source paths, presence, index status, and bytes; it does not claim an unavailable project-wide SNV/indel VCF exists | Section 7.6; `results/data_delivery_availability_manifest.tsv`; `results/customer_request_coverage.tsv` |

### 8.3 Limitations

This work is a deeper analysis of the same short-read WGS alignments; it is not an independent biological replicate or assay. A 2-bp overlap in split-read alignment representation is consistent with microhomology but cannot by itself prove a unique breakpoint in a repetitive locus. Manta/DELLY agreement, if it had been observed, would also have been algorithmic evidence from the same data rather than an external confirmation.

No independent breakpoint-spanning polymerase chain reaction (PCR) result, new sequencing dataset, or raw fluorescence in situ hybridization (FISH) image/scoring table was supplied. The available local depth is below the 60x sensitivity context used for the 3-fragment/5% criterion. Therefore, zero support in the parental/control samples is not a validated high-sensitivity negative result. No provenance-qualified project-wide SNV/indel VCF source was identified in the current inventory.

### 8.4 Next steps

The scientific next evidence for C2 is to test the previously selected in-silico breakpoint-spanning PCR pair on independent DLD1 genomic DNA with appropriate controls and measured analytical sensitivity. If a separate sensitivity study is requested, its measured local short-read WGS depth should be evaluated against the pre-specified 60x context before interpreting a three-fragment/5% detection target. Raw FISH images and a per-cell scoring table would be needed before linking a genomic candidate to a cytogenetic subpopulation.

No external data transfer, client contact, commercial quote, or billing action was initiated in preparing this scientific report. The client’s stated $840 cap and any transfer method require separate Athenomics administrative confirmation before an external action or a scope beyond the analysis reported here.

## 9. Deliverable Files

| Path in this supplementary package | Contents | Purpose |
| :--- | :---: | :--- |
| `{REPORT_NAME}` | This English report | Entry point for scientific review |
| `figures/{FIGURE_STEM}.png` and `.pdf` | One client-question completion figure in two formats | Visual summary of the P41-to-P42 refinement, rescue impact, and caller comparison |
| `results/c2_junction_resolution_summary.tsv` | One C2 summary row | Machine-readable microhomology/gap summary |
| `results/c2_read_end_junction_resolution.tsv` | Five resolved read-end rows | Alignment-level evidence retained with FLAG, CIGAR, MAPQ, and query coordinates |
| `results/candidate_manta_delly_concordance.tsv` | Twenty candidate-by-DLD1-sample comparisons | Exact two-end and anchor-only caller classifications |
| `results/clone_specific_comparison_with_criteria.tsv` | Twenty candidate-by-DLD1-sample rows | Support, depth, and gate criteria |
| `results/read_rescue_candidate_impact.tsv` | Four C1-C4 rows | Before/after rescue accounting |
| `results/customer_request_coverage.tsv` | Five request-to-output rows | Scope coverage and interpretation boundaries |
| `results/data_delivery_availability_manifest.tsv` | 109 source-inventory rows | Transfer-planning inventory without duplicating source data |
| `results/customer_question_summary.tsv` and `results/claim_evidence_matrix.tsv` | Compact answer and claim-source tables | Review-friendly entry points to the detailed evidence |
| `methods/source_manifest.json` | Source SHA-256 and inventory summary | Reproducibility and source traceability |

---

Prepared for scientific review; not sent to the client.

Zhen Gao, PhD  
Principal Bioinformatics Scientist  
Athenomics
"""


def main() -> None:
    args = parse_args()
    out_dir = args.out_dir.resolve()
    delivery = args.delivery.resolve()
    if delivery.exists():
        raise FileExistsError(f"Refusing to overwrite existing delivery directory: {delivery}")
    required_names = [
        "c2_junction_resolution_summary.tsv",
        "c2_read_end_junction_resolution.tsv",
        "candidate_manta_delly_concordance.tsv",
        "clone_specific_comparison_with_criteria.tsv",
        "read_rescue_candidate_impact.tsv",
        "customer_request_coverage.tsv",
        "data_delivery_availability_manifest.tsv",
        "customer_request_audit.json",
    ]
    figures = out_dir / "figures"
    required = [out_dir / name for name in required_names] + [
        figures / f"{FIGURE_STEM}.png", figures / f"{FIGURE_STEM}.pdf"
    ]
    missing = [str(path) for path in required if not path.is_file()]
    if missing:
        raise FileNotFoundError("Missing required new analysis output(s):\n" + "\n".join(missing))

    delivery.mkdir()
    results = delivery / "results"
    figure_dir = delivery / "figures"
    methods = delivery / "methods"
    results.mkdir()
    figure_dir.mkdir()
    methods.mkdir()

    copies: dict[Path, Path] = {
        **{out_dir / name: results / name for name in required_names},
        figures / f"{FIGURE_STEM}.png": figure_dir / f"{FIGURE_STEM}.png",
        figures / f"{FIGURE_STEM}.pdf": figure_dir / f"{FIGURE_STEM}.pdf",
    }
    for source, destination in copies.items():
        shutil.copy2(source, destination)

    c2 = read_tsv(out_dir / "c2_junction_resolution_summary.tsv")[0]
    rescue = read_tsv(out_dir / "read_rescue_candidate_impact.tsv")
    caller = read_tsv(out_dir / "candidate_manta_delly_concordance.tsv")
    clone = read_tsv(out_dir / "clone_specific_comparison_with_criteria.tsv")
    availability = read_tsv(out_dir / "data_delivery_availability_manifest.tsv")
    manifest = source_manifest(availability)
    if (c2["independent_qname_count"], c2["split_read_end_count"], c2["query_overlap_bp_distribution"],
            c2["unmapped_query_gap_bp_distribution"]) != ("3", "5", "2bp:5", "0;0;0;0;0"):
        raise ValueError("C2 summary differs from the approved customer-question analysis result")
    if len(caller) != 20 or sum(
        int(int(row["manta_exact_two_end_record_count"]) > 0 and int(row["delly_exact_two_end_record_count"]) > 0)
        for row in caller
    ) != 0:
        raise ValueError("Caller-concordance source is incomplete or inconsistent")
    if len(rescue) != 4 or any(row["raw_dual_endpoint_qname_count_before_rescue"] != row["raw_dual_endpoint_qname_count_after_rescue"] for row in rescue):
        raise ValueError("Rescue-impact source is incomplete or inconsistent")
    if len(clone) != 20:
        raise ValueError("DLD1 clone-comparison source is incomplete")

    summary_rows = [
        {"topic": "P41-to-P42 C2 refinement", "answer": "P41 EV0623 crosswalks to C2; 3 QNAMEs and 5 split-read ends carry 2-bp AT overlap with 0-bp query gap.", "boundary": "Repeat-context representation and same-WGS evidence are not independent biological confirmation."},
        {"topic": "Manta/DELLY exact two-end concordance", "answer": "0/20 Manta, 0/20 DELLY, and 0/20 both-caller exact two-end matches.", "boundary": "Five C4 rows are DELLY anchor-only/nonconcordant records, not BND support."},
        {"topic": "Read-rescue impact", "answer": "Five unmapped-mate QNAMEs were retained; C1-C4 each had zero overlap and unchanged counts.", "boundary": "Zero impact is candidate-specific and does not discard the rescued evidence."},
        {"topic": "DLD1 clone comparison", "answer": "No C1-C4 candidate meets the high-confidence clone/doxycycline gate.", "boundary": "No validated allele fraction is reported; available control depth is below 60x sensitivity context."},
        {"topic": "Source-file availability", "answer": f"{manifest['BAM']['existing_files']}/18 BAMs and indexes are present; source categories are inventoried.", "boundary": "No provenance-qualified project-wide SNV/indel VCF source was identified."},
    ]
    write_tsv(results / "customer_question_summary.tsv", summary_rows, ["topic", "answer", "boundary"])
    claim_rows = [
        {"claim": "C2 read-end representation is consistent with 2-bp AT microhomology and no observed query gap", "source_file": "c2_read_end_junction_resolution.tsv; c2_junction_resolution_summary.tsv", "interpretation_boundary": "Not independent confirmation or a unique repeat-context breakpoint proof"},
        {"claim": "No exact two-end Manta/DELLY match for C1-C4 in the DLD1 comparison", "source_file": "candidate_manta_delly_concordance.tsv", "interpretation_boundary": "Caller comparison is algorithmic evidence from the same WGS dataset"},
        {"claim": "C1-C4 rescue counts are unchanged", "source_file": "read_rescue_candidate_impact.tsv", "interpretation_boundary": "Only addresses overlap of five rescued unmapped mates with C1-C4"},
        {"claim": "No C1-C4 high-confidence clone/doxycycline-specific event", "source_file": "clone_specific_comparison_with_criteria.tsv", "interpretation_boundary": "Requires all support/depth/sequence gates; zero in low-depth controls is not exclusion"},
        {"claim": "BAM/VCF/BCF sources have an availability inventory; SNV/indel VCF is not claimed", "source_file": "data_delivery_availability_manifest.tsv", "interpretation_boundary": "Inventory only; no external transfer was initiated"},
    ]
    write_tsv(results / "claim_evidence_matrix.tsv", claim_rows, ["claim", "source_file", "interpretation_boundary"])

    report = make_report(manifest)
    (delivery / REPORT_NAME).write_text(report)
    provenance = {
        "analysis_date": date.today().isoformat(),
        "scope": "Incremental client-question completion for P41 Section 7.3; no P41/P42 delivery was modified or duplicated.",
        "analysis_source": str(out_dir),
        "delivery": str(delivery),
        "copied_source_sha256": {str(source.relative_to(out_dir)): sha256(source) for source in copies},
        "inventory_summary": manifest,
        "report": REPORT_NAME,
    }
    (methods / "source_manifest.json").write_text(json.dumps(provenance, indent=2) + "\n")


if __name__ == "__main__":
    main()
