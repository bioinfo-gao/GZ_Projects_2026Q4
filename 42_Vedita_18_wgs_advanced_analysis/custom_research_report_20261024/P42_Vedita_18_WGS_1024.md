# P42 whole-genome sequencing advanced analysis and P41 Section 7.3 follow-up

**Report Date:** 2026-10-24  
**Prepared by:** Zhen Gao, PhD, Principal Bioinformatics Scientist, Athenomics  
**Analysis Platform:** Linux high-performance computing server  
**Reference Genome:** *Homo sapiens* GRCh38  
**Sample context:** DLD1 parental population and clone samples, with and without doxycycline, as named in the follow-up request; tissue of origin was not specified by the client.  
**Sequencing data context:** The WGS data amount is approximately 25-30x for most samples. The restricted primary/non-duplicate/MAPQ>=20 depth used in this candidate comparison is reported separately as an analysis denominator, not as the WGS service-yield metric.  
**Delivery status:** Prepared for scientific review and client delivery; not sent to the client.

## 1. Objectives

This single delivery consolidates the complete P42 analysis and the follow-up research requested above P41 report Section 7.3. It provides the full 18-sample short-read whole-genome sequencing (WGS) evidence review, refined C2 junction analysis, caller concordance, completed unmapped-mate-rescue impact, DLD1 control-versus-doxycycline comparison, an in-silico C2 assay design, and the requested transferable source files. The comparison framework explicitly separates service-level WGS data yield from the filtered standard-evidence depth used to interpret candidate observations.

## 2. Core conclusions

1. P42 evaluated 18 WGS samples across the requested chr11 repeat-region evidence window. The completed inventory retained low-MAPQ, multi-mapping, secondary, supplementary, soft-clipped, and unmapped-mate contexts rather than treating their absence from a filtered candidate list as biological absence.
2. P41 EV0623 is carried forward as P42 C2, not a new event: P41 reported chr11:114,607,192 to chr21:20,648,241, whereas the P42 representation is chr11:114,607,191 to chr21:20,648,242. The one-base difference at each end is a representation difference for the same geometry.
3. C2 is now resolved to five split-read ends from three distinct QNAME labels. All five have a 2-base-pair (bp) `AT` query overlap and a 0-bp query gap at MAPQ 60. This is internally consistent alignment evidence, not independent biological confirmation in a repeat context.
4. Across 20 fixed candidate-by-DLD1-sample comparisons, exact two-end concordance is 0/20 for Manta and 0/20 for DELLY. Five C4 DELLY rows are chr11-only anchors without the chr3 partner and are not called concordant.
5. C1-C4 remain candidate observations; none establishes a clone- or doxycycline-specific rearrangement. Their raw and standard-eligible QNAME counts are reported separately, without treating a support/depth ratio as an allele fraction or a comparator zero as an absence call. No retrospective depth pass/fail gate is used.
6. The C2 assay pair is an in-silico design only. The requested 18 BAM files and indexes, Manta VCFs, and DELLY records are included here. A provenance-qualified project-wide SNV/indel VCF was not identified and is explicitly unavailable.

## 3. Study inputs and sample information

The analysis includes 18 named WGS samples. The WGS data amount is approximately 25-30x for most samples. For the standard candidate-comparison denominator only, base-weighted autosomal depth from primary, non-duplicate, non-QC-fail, non-supplementary alignments with MAPQ>=20 ranges from 8.513417x to 20.125928x across the cohort. This restricted denominator is not a measurement of WGS service yield. The DLD1 follow-up comparison uses S01-S05.

| Sample | Requested role | Treatment label | Standard-evidence autosomal depth | Use in P41 Section 7.3 follow-up |
| :--- | :---: | :---: | :---: | :--- |
| S01_DLD1_pop | Parental population | Not stated as +doxycycline | 15.449764x | Comparison reference |
| S02_DLD1_c1_ctrl | Control clone | Control | 13.677017x | Comparison reference |
| S03_DLD1_c2_ctrl | Control clone | Control | 15.172271x | Comparison reference |
| S04_DLD1_c3_dox | DLD1 clone | +Doxycycline | 15.415246x | Case group |
| S05_DLD1_c4_dox | DLD1 clone | +Doxycycline | 17.823929x | Case group |

All 18 BAM files and their indexes are included in `genome_wide_data/bam/`. The source BAM total is 907,440,957,244 bytes (907.44 GB; 845.12 GiB). Client BAM headers have local execution paths removed; alignment records and their per-file counts are retained and audited.

## 4. Analysis design and decision criteria

| Research question | Evidence unit and retained context | Decision criterion | Permitted conclusion |
| :--- | :---: | :---: | :--- |
| What depth is used in this candidate comparison? | Base-weighted chr1-chr22 depth from primary, non-duplicate, non-QC-fail, non-supplementary alignments with MAPQ>=20 | Report this restricted denominator separately from the approximately 25-30x WGS data-yield context | Analysis-specific evidence context; not a service-yield or pass/fail metric |
| Does a window contain mapping evidence? | Original alignment records, retaining primary, secondary, supplementary, soft-clipped, multi-mapping, and unmapped-mate context | Evidence counts and filter-impact audit before candidate restriction | Mapping evidence, not an independent molecular count |
| Is C2 read-end structure internally consistent? | Sequence-bearing records grouped by QNAME and read end, retaining original FLAG, CIGAR, MAPQ, and coordinates | Query overlap/gap and overlap sequence across resolved ends | Consistent split-read representation, not a uniquely proven breakpoint |
| Do Manta and DELLY match C1-C4? | All audited records from 18 samples and four caller stages | Both BND ends within 50 bp in either order | Algorithmic concordance from the same WGS data only |
| Is clone/doxycycline specificity established? | Fixed-candidate raw and standard-eligible QNAME counts, local depth context, sequence-level evidence, and independent evidence if supplied | >=3 standard-eligible QNAMEs is a support-strength descriptor for prioritizing a positive candidate review, not an absence gate. A positive repeat-context conclusion also needs coherent target evidence and independent confirmation. | Candidate status or an established positive conclusion; a comparator zero at the available depth is reported as not observed, not absent |
| Did rescue alter C1-C4? | Recovered unmapped-mate QNAMEs intersected with candidate QNAMEs | Before/after candidate counts and overlap | Measured rescue impact on C1-C4 |

No direct allele-fraction cutoff is applied: raw support divided by the restricted standard-evidence depth is not a validated allele-fraction estimator in this repeat context. The client-facing depth tables make the denominator explicit and retain the boundary on any absence-style interpretation.

## 5. Methods

P42 evaluated 59 assemblies yielding 90 contigs and 32 alternative placements. From 282 candidate observations, 94 were parseable under the fixed evidence model; four case-only candidates (C1-C4) were retained for detailed review. The raw evidence ledger preserves alignment context and records filter impact rather than allowing a prior caller table to define the data universe.

For the P41 Section 7.3 follow-up, P41 EV0623 was crosswalked to C2. Sequence-bearing alignment records were resolved at `(QNAME, read end)` level using CIGAR-derived query intervals. Manta and DELLY records were compared with each fixed candidate under a two-end, 50-bp tolerance rule. A single nearby anchor, a different partner, or a different SV type remained a retained nonconcordant record. The rescue analysis intersected the completed unmapped-mate ledger with the C1-C4 QNAME sets. The C2 assay is an in-silico primer design and specificity screen; no new assay was performed.

For this revision, the standard-evidence autosomal denominator was recalculated from the existing 100-kb bins on chr1-chr22, using the same primary/non-duplicate/non-QC-fail/non-supplementary/MAPQ>=20 track for every sample. It is used solely to qualify interpretation of an observed zero in the candidate comparison. It is not substituted for the WGS data-yield context and is not used as a service-quality conclusion.

## 6. Results

### 6.1 Standard-evidence depth context

The WGS data amount is approximately 25-30x for most samples. The restricted standard-evidence autosomal denominator used in this comparison ranges from 8.513417x to 20.125928x over all 18 samples and from 13.677017x to 17.823929x in DLD1 S01-S05. These two metrics are intentionally not interchangeable. `results/standard_evidence_autosomal_depth.tsv` and `results/standard_evidence_depth_summary.tsv` give the calculation and definition.

### 6.2 Complete P42 evidence inventory

The full evidence inventory and assemblies support a repeat-aware interpretation. C1-C4 remain candidate mapping observations rather than established clone- or doxycycline-specific structural variants. Counts, placements, and candidate-level support are supplied in `results/` and Figures 1-12.

### 6.3 P41 Section 7.3 refinement: C2 read-end junction structure

P41 EV0623 maps to the P42 C2 geometry. Three QNAME labels yielded five resolved split-read ends. Every end has MAPQ 60 for both mapped blocks, a 2-bp `AT` query overlap, and no unaligned query gap. This narrows the technical description beyond the P41 split-read count while preserving the limitation that a repeat-context short-read alignment is not an independent biological confirmation.

### 6.4 Manta and DELLY comparison

| Candidate | Exact Manta two-end matches across S01-S05 | Exact DELLY two-end matches across S01-S05 | Both callers exact | Anchor-only/nonconcordant rows |
| :--- | :---: | :---: | :---: | :---: |
| C1 | 0 | 0 | 0 | 0 |
| C2 | 0 | 0 | 0 | 0 |
| C3 | 0 | 0 | 0 | 0 |
| C4 | 0 | 0 | 0 | 5 |
| Total | 0/20 | 0/20 | 0/20 | 5/20 |

The five C4 nonconcordant rows are recurrent DELLY chr11 insertion anchors without the chr3 partner. They are retained as caller context and are not promoted to evidence for a chr11-chr3 BND.

### 6.5 Unmapped-mate rescue impact

Five unmapped-mate QNAMEs were recovered across the full locus. None overlaps a C1-C4 dual-endpoint QNAME; raw support is unchanged: C1 1 to 1, C2 3 to 3, C3 2 to 2, and C4 1 to 1.

### 6.6 DLD1 control versus doxycycline comparison

| Candidate | Raw support in S04 / S05 | Standard-eligible support in S04 / S05 | Raw support in S01-S03 | Interpretation |
| :--- | :---: | :---: | :---: | :--- |
| C1 | 0 / 1 | 0 / 1 | 0 / 0 / 0 | Low-support S05 candidate; clone/doxycycline specificity is not established |
| C2 | 0 / 3 | 0 / 0 | 0 / 0 / 0 | P41 carry-forward candidate; clone/doxycycline specificity is not established |
| C3 | 1 / 1 | 0 / 1 | 0 / 0 / 0 | Low-support observations in both +doxycycline samples; specificity is not established |
| C4 | 0 / 1 | 0 / 1 | 0 / 0 / 0 | Low-support S05 candidate; clone/doxycycline specificity is not established |

For C2, S05 has a raw-support/depth proxy of 3 / 14.760 = 0.203252. It is not a validated allele fraction. C2 standard local depths in S01-S03 are 21.146x, 19.126x, and 23.101x. The zero raw counts in those comparators are therefore reported as not observed within the available standard-evidence context; they are not used as a high-sensitivity exclusion or as evidence against clone/doxycycline specificity.

### 6.7 In-silico C2 assay design

The selected C2 junction-spanning primer pair is `GTAATCACGGTTACCCAGTGGAT` and `AGCAATAGTTCTCAACAGTGGGT`, with an expected 274-bp product on the virtual junction template. This is an in-silico design and specificity screen only; it does not establish assay performance or confirm C2.

## 7. Figures and supplied results

Figures 1-12 document the complete P42 evidence, assembly, placement, and candidate review. Figure 13 documents the C2 in-silico assay design. Figure 14 summarizes the P41 Section 7.3 client-question follow-up. Every figure is provided in PDF and PNG in `figures/`; `figures/figure_manifest.json` records the 14 formal figures. Machine-readable result tables are in `results/`, including the evidence crosswalk, C2 read-end ledger, caller comparison, rescue impact, clone comparison, standard-evidence depth context, assay selection, and client-request coverage.

## 8. Conclusions, Limitations and Next Steps

### 8.1 What this analysis establishes

The evidence establishes a precise technical characterization of C2, the absence of exact two-end caller concordance under the stated rule, the measured zero C1-C4 rescue increment, and the fact that C1-C4 do not establish clone- or doxycycline-specific rearrangements. It also clearly separates the WGS data-yield context from the restricted depth denominator used to interpret standard candidate evidence, and supplies the specified source-data file categories in one package.

### 8.2 Practical reading: claim-to-evidence map

| Claim | Direct evidence | Interpretation boundary |
| :--- | :---: | :--- |
| C2 refines P41 EV0623 | `results/c2_junction_resolution_summary.tsv`; `results/c2_read_end_junction_resolution.tsv` | Same WGS alignment evidence, not independent validation |
| C2 has 2-bp `AT` overlap and 0-bp gap | Five resolved read ends in `results/c2_read_end_junction_resolution.tsv` | Repeat-context placement remains ambiguous |
| No exact Manta or DELLY support | `results/candidate_manta_delly_concordance.tsv` | Caller concordance would still be algorithmic evidence |
| Rescue did not increase C1-C4 support | `results/read_rescue_candidate_impact.tsv` | Does not exclude rescue value elsewhere |
| Clone/doxycycline specificity is not established for C1-C4 | `results/clone_specific_comparison_with_criteria.tsv`; `results/standard_evidence_depth_summary.tsv` | Comparator zeros are not absence calls; no retrospective depth gate or direct allele-fraction threshold is applied |
| C2 assay is design-only | `results/assay_design_selection.tsv`; `methods/virtual_junction_templates.fasta` | No independent assay result is included |

### 8.3 Limitations

This is a deeper analysis of the same short-read WGS data, not an independent biological replicate. A split-read overlap supports a mapping representation but cannot alone prove a unique breakpoint in a repetitive locus. The filtered standard-evidence depth is intentionally distinct from WGS data yield and constrains only interpretation of a comparator zero. No independent breakpoint-spanning assay, new sequencing, raw fluorescence in situ hybridization input, or scored cytogenetic result was supplied for this package. The unavailable project-wide SNV/indel VCF is not substituted by Manta candidate small-indel output.

### 8.4 Next evidence

The appropriate next evidentiary step for C2 is experimental assessment of the supplied junction-spanning assay design and, if available, integration of independently generated cytogenetic scoring with pre-specified interpretation rules. Any later low-fraction sensitivity study should define its sample type, assay recovery, and interpretation rule prospectively rather than retrospectively imposing a different WGS depth gate. Any later analysis should retain the raw-alignment and filter-impact audit used here.

## 9. Delivery contents

| Directory or file | Contents | Status |
| :--- | :---: | :--- |
| `P42_Vedita_18_WGS_1024.md` | This integrated client-facing report | Included |
| `figures/` | 14 formal figures, each as PDF and PNG | Included |
| `results/` | Machine-readable analysis and follow-up tables | Included |
| `results/standard_evidence_autosomal_depth.tsv` | All-sample restricted autosomal analysis denominator, separately defined from WGS data yield | Included |
| `results/standard_evidence_depth_summary.tsv` | Cohort and DLD1 summary of the same restricted denominator | Included |
| `methods/` | Client-safe methods manifest, `depth_interpretation_method.tsv`, and virtual junction template | Included |
| `genome_wide_data/bam/` | 18 BAM files and 18 BAM indexes; alignment records preserved and local execution paths removed from headers | Included |
| `genome_wide_data/manta/` | Filtered SV, candidate SV, and candidate small-indel VCFs; variant records preserved and local execution paths removed from headers | Included |
| `genome_wide_data/delly/` | Discovery and region-genotyped VCFs, converted from source BCF with indexes; variant records preserved and local execution paths removed from headers | Included |
| `results/data_inventory.tsv` | Client-relative inventory of all supplied data | Included |
| `methods/vcf_header_sanitization_audit.tsv` | Per-file variant-body hash comparison and removed-header-path count | Included |
| `methods/bam_header_sanitization_audit.tsv` | Per-file alignment-count comparison and removed-header-path count | Included |
| Project-wide SNV/indel VCF | No provenance-qualified source identified | Unavailable; not represented as present |
