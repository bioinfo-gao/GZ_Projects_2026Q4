#!/usr/bin/env python3
"""Write the one-source client report from the P42 rescue result tables."""
import csv, pathlib, sys

samples_csv, results_tsv, summary_tsv, output = map(pathlib.Path, sys.argv[1:5])
samples = list(csv.DictReader(samples_csv.open()))
rows = list(csv.DictReader(results_tsv.open(), delimiter="\t"))
summary = list(csv.DictReader(summary_tsv.open(), delimiter="\t"))
assert len(rows) == 179 and len(samples) == 10 and len(summary) == 10
status = {x: sum(r["rescue_status"] == x for r in rows) for x in
          ("multi_fragment_supported", "single_fragment_supported", "ambiguous_reconstruction", "not_reconstructed", "not_evaluable")}
unreconstructed = [r for r in rows if r["rescue_status"] == "not_reconstructed"]
assert sum(status.values()) == 179 and len(unreconstructed) == 2
fmt = lambda n: f"{n:,}"
unrec_text = "; ".join(f"{r['candidate_id']} ({r['type']}, {r['chr1']}:{r['pos1']} to {r['chr2']}:{r['pos2']})" for r in unreconstructed)
sample_rows = "\n".join(f"| {r['sample']} | {r['candidates']} | {r['multi_fragment_supported']} | {r['single_fragment_supported']} | {r['ambiguous_reconstruction']} | {r['not_reconstructed']} |" for r in summary)
input_rows = "\n".join(f"| {s['sample_id']} | {s['client_label']} | {s['cohort_role']} | {s['source_candidate_count']} |" for s in samples)
text = f'''# P42 Supplementary chr11 Junction Read-Rescue Report

**Report Date:** @@REPORT_DATE@@  
**Prepared by:** Zhen Gao, PhD, Principal Bioinformatics Scientist, Athenomics  
**Analysis Platform:** Linux high-performance computing server  
**Species:** not specified by client; GRCh38 was used as the analytical reference coordinate system  
**Tissue / Cell:** not specified by client

## 1. Objectives

This supplementary analysis extends the previously delivered chr11 report with read-level rescue for the client-defined target interval, chr11:114,604,000-114,625,500. It re-extracted the sequences from the original candidate-supporting fragments and independently re-aligned them to the same GRCh38 reference with all alternative alignments retained. The scope contains S01-S10, including S02, S03 and S10; Mutu I and DG75 samples are not part of this supplementary analysis.

## 2. Core Conclusions

1. **The two recorded breakpoint sides were independently reconstructed for {fmt(status['multi_fragment_supported'] + status['single_fragment_supported'])} of {fmt(len(rows))} candidate junctions.** This comprises {status['multi_fragment_supported']} candidates supported by two or more uniquely re-aligned fragments and {status['single_fragment_supported']} supported by one uniquely re-aligned fragment.
2. **Most reconstructed candidates remain single-fragment observations.** A reconstructed alignment shows that the original read sequence can again form the recorded junction; it does not by itself validate an underlying structural alteration.
3. **Two candidates were not reconstructed by their extracted sequence on independent re-alignment:** {unrec_text}. They should not be treated as read-rescue-supported.
4. **No candidate is called confirmed solely from this supplementary analysis.** Events requiring a definitive sequence-level conclusion need targeted manual breakpoint curation or an orthogonal confirmation assay; this contracted read-rescue scope did not include manual local assembly or breakpoint adjudication.

## 3. Key Findings

The pre-registered candidate definition required at least one non-duplicate P41 split read with mapping quality (MAPQ) >=20 on both P41 segments. This yielded {fmt(len(rows))} candidates across the ten requested samples. Independent `bwa-mem2 -a` re-alignment reconstructed both recorded breakpoint sides within 20 bp for {fmt(status['multi_fragment_supported'] + status['single_fragment_supported'])} candidates. No reconstructed candidate was classified as ambiguous under the defined MAPQ criterion; {status['not_reconstructed']} were not reconstructed and {status['not_evaluable']} lacked an evaluable source read.

| Sample | Candidates | Multi-fragment supported | Single-fragment supported | Ambiguous reconstruction | Not reconstructed |
| :--- | :---: | :---: | :---: | :---: | :---: |
{sample_rows}

## 4. Sample Information

Sample labels and scope were taken from the current client correspondence. The cohort labels below identify the requested comparison context and are not additional sample metadata.

| Sample | Client label | Cohort role | Candidates |
| :--- | :---: | :---: | :---: |
{input_rows}

## 5. Analysis Rationale and Decision Criteria

P41 had reported candidate junctions in a repeat-affected region. Read rescue therefore tested whether the original supporting read sequence can independently reconstruct the same two breakpoint ends while exposing alternative alignments. A candidate was retained before re-alignment only if P41 had at least one non-duplicate split read with MAPQ >=20 on both segments. For each retained source fragment, an independently re-aligned fragment was called **uniquely reconstructed** only when one alignment reached each recorded breakpoint end within 20 bp and both matching alignments had MAPQ >=20. Candidate-level labels were: `multi_fragment_supported` (at least two uniquely reconstructed fragments), `single_fragment_supported` (one), `ambiguous_reconstruction` (reconstructed but failing the MAPQ condition), `not_reconstructed`, or `not_evaluable`.

## 6. Methods

The coordinate-sorted, duplicate-marked Binary Alignment/Map (BAM) file for each in-scope sample was queried only within chr11:114,554,000-114,675,500; the source BAM files were not modified. Candidate support-fragment sequences were re-extracted and independently aligned with `bwa-mem2` 2.2.1 using `mem -a -Y` against the Genome Analysis Toolkit (GATK) GRCh38 analysis-set reference. The `-a` option retains alternative alignments, preventing a repeat-derived read from being presented as uniquely mapped by omission. The resulting alignments were evaluated against each recorded P41 breakpoint coordinate using the 20 bp endpoint tolerance described above.

## 7. Results

All {fmt(len(rows))} retained candidates had at least one source fragment re-extracted from the corresponding BAM. The distribution was {status['multi_fragment_supported']} multi-fragment supported, {status['single_fragment_supported']} single-fragment supported, {status['ambiguous_reconstruction']} ambiguous reconstruction, {status['not_reconstructed']} not reconstructed, and {status['not_evaluable']} not evaluable. The complete per-candidate coordinates, original P41 support counts, rescued-fragment counts, and classifications are provided in `junction_rescue/rescue_candidate_results.tsv`.

The evidence strength differs materially between the two supported categories: {status['multi_fragment_supported']} candidates have at least two independently re-aligned support fragments, whereas {status['single_fragment_supported']} have one. The latter remain observations from a single sequenced molecule and should not be interpreted as validated structural variants without the additional confirmation described below.

## 8. Conclusions, Limitations and Next Steps

### 8.1 Practical conclusions by evidence category

| Evidence category | Number | Practical reading |
| :--- | :---: | :--- |
| Multi-fragment supported | {status['multi_fragment_supported']} | Stronger read-level support; select individual candidates for manual curation only if they are biologically relevant. |
| Single-fragment supported | {status['single_fragment_supported']} | Sequence-level reconstruction is reproducible, but one molecule is insufficient for structural-variant confirmation. |
| Not reconstructed | {status['not_reconstructed']} | Do not use as read-rescue-supported evidence. |

### 8.2 Practical reading with evidence sources

The classification table `junction_rescue/rescue_candidate_results.tsv` is the evidence source for all counts in Section 8.1. `multi_fragment_supported` requires at least two uniquely reconstructed fragments under the Section 5 criteria. `single_fragment_supported` requires one. The two non-reconstructed candidates are identified in Section 2 and remain in the table to preserve the complete candidate record rather than being silently removed.

### 8.3 Limitations

This analysis independently re-evaluates the same short-read sequences; it is not independent biological replication. It does not perform local assembly, base-pair breakpoint curation, or manual visual review. In a repeat-affected region, a single uniquely re-aligned short fragment is useful read-level evidence but does not determine clonality, allele fraction, or the full structural configuration.

### 8.4 Suggested confirmation

For any candidate selected for a biological decision, review its local alignments and use a targeted orthogonal assay designed across the proposed junction. Prioritise the four multi-fragment-supported candidates before single-fragment observations.

## 9. Deliverable Files

| File | Description |
| :--- | :--- |
| `junction_rescue/rescue_candidate_results.tsv` | All {fmt(len(rows))} candidates with P41 coordinates, original support, re-extracted fragment count, rescue counts and final read-rescue label. |
| `junction_rescue/rescue_sample_summary.tsv` | Per-sample totals and read-rescue categories for S01-S10. |
| `junction_rescue/candidate_counts.tsv` | Candidate-count audit trail for the ten-sample scope. |

Zhen Gao, PhD  
Principal Bioinformatics Scientist  
Athenomics
'''
output.write_text(text)
print(f"wrote report source for {len(rows)} candidates")
