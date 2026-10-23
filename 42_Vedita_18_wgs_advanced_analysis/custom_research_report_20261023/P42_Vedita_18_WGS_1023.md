# P42 Full-Cohort Whole-Genome Sequencing Analysis Report

**Report Date:** 2026-10-23  
**Prepared by:** Zhen Gao, PhD, Principal Bioinformatics Scientist, Athenomics  
**Analysis Platform:** Linux high-performance computing server  
**Species:** *Homo sapiens*; Genome Reference Consortium Human Build 38 (GRCh38) was the analytical reference coordinate system  
**Tissue / Cell:** not specified by client

## 1. Objectives

This full-cohort whole-genome sequencing (WGS) analysis examined all 18 supplied alignment datasets (S01–S18). It assessed chromosome 11 (chr11) coverage, raw junction evidence, local assemblies, alternative genomic placements, and external structural-variant callers. The primary comparison used the analysis labels S04/S05 (`+dox`, doxycycline treatment) against S01 (parental) and S02/S03 (control clones). The remaining 13 samples were retained as cohort background rather than treated as DLD1 biological controls.

The analysis was designed to determine whether the reported chr11/FISH observation could be resolved to a Dox- or clone-specific structural breakpoint. It extends the earlier supplementary read-rescue package with an all-sample, repeat-aware evidence model; the prior delivery remains an intact historical record.

## 2. Core Conclusions

1. **A high-confidence Dox- or clone-specific structural breakpoint was not established.** Four caller-defined patterns had raw two-endpoint evidence only in the `+dox` samples, but none met the pre-specified confirmation gate.
2. **The four patterns should be regarded as candidates, not confirmed events.** Their control samples had 14.7–23.1× standard local depth, below the 60× depth required for a three-fragment / 5% sensitivity context; the available assemblies do not demonstrate a sequence-level junction.
3. **The chr11 coverage pattern is a useful screening observation but does not identify a breakpoint.** S04/S05 show lower normalized depth across substantial parts of the discovery window, but coverage alone cannot distinguish copy-number change from repeat mappability or library-depth effects.
4. **The analysis preserves rather than suppresses repeat-associated ambiguity.** Alternative placements were observed for 32 of 90 local-assembly contigs. This prevents a primary alignment label from being over-interpreted as a unique genomic location.

Section 3 presents the measured findings by theme. Section 8.2 provides the evidence files and decision basis behind each practical conclusion; Section 8.3 states what the data cannot determine, and Section 8.4 identifies the focused confirmation work for the four candidates.

## 3. Key Findings

1. **Cohort coverage and raw evidence were available for all 18 samples.** Canonical genome mean depth ranged from 8.2× to 19.3×. The chr11 discovery-window inventory retained primary, secondary, supplementary and mapping-quality-zero records as separate auditable categories.
2. **Local assembly completed for all 59 selected loci.** The assemblies yielded 90 contigs; 32 contigs had more than one genomic placement. Contig mapping is therefore context for breakpoint assessment, not confirmation by itself.
3. **External callers produced a broad candidate universe, not a final event list.** Genome Rearrangement IDentification Software Suite (GRIDSS) reported 278 normalized target-related breakend events across the 18-sample caller run. Structural Variation Analysis by Assembly (SvABA) reported four germline breakend events and no somatic breakend events in the corrected five-sample DLD1 run.
4. **Fixed-candidate raw re-genotyping separates geometry from evidence strength.** Of 282 caller breakend events, 94 had a uniquely parsed two-endpoint geometry. Seven showed raw two-endpoint support in at least one case and control, 21 only in controls, and four only in `+dox` samples. None of the four case-only patterns passes the coverage and sequence-level assembly requirements for a high-confidence conclusion.

## 4. Sample Information

The client did not provide an independent tissue/cell metadata field with this delivery. The table below records the supplied analysis labels and comparison roles; labels are not used to infer additional biological metadata.

| Analysis samples | Supplied analysis labels | Role in this report |
| :--- | :---: | :---: |
| S01–S05 | DLD1 population / control / `+dox` labels | Primary fixed-candidate comparison: S04/S05 versus S01/S02/S03 |
| S06–S10 | GM12878 population / clone labels | Cohort background and repeat-context comparison |
| S11–S14 | Mutu I population / clone labels | Cohort background and repeat-context comparison |
| S15–S18 | DG75 population / clone labels | Cohort background and repeat-context comparison |

## 5. Analysis Rationale and Decision Criteria

The chr11 discovery interval was chr11:114,604,000–114,625,500 (GRCh38, 1-based coordinates), expanded by 50 kilobases on each side for evidence discovery. In this repeat-affected setting, mapping records were not discarded merely because they were secondary, supplementary, low mapping quality, or duplicate-marked. Instead, their alignment class and mapping context were retained.

An independent-fragment count was defined as one unique query template name (QNAME) within one sample and fixed candidate. Multiple records or read ends from that QNAME remain in the audit ledger but are not counted as separate fragments. A raw two-endpoint fragment required either a supplementary-alignment (SA) tag or discordant paired-end mate to match both caller endpoints within 50 base pairs (bp). Soft-clipped or supplementary records matching only one endpoint were recorded as anchor context, not as complete breakpoint support.

The pre-specified high-confidence gate required at least three standard-eligible independent fragments (primary, non-duplicate, non-quality-control-fail, mapping quality at least 20) and a 5% sensitivity context. At the three-fragment threshold, 5% requires at least 60× local depth. A sequence-level contig spanning the junction would also be required; a contig that only overlaps one chr11 anchor is not sufficient.

## 6. Methods

Coordinate-sorted, duplicate-marked Binary Alignment/Map (BAM) files were queried read-only against the same GRCh38 reference. The full chr11 discovery window was inventoried before fixed candidates were considered. Soft clips, supplementary alignments, discordant mates and unmapped-mate rescue were retained in machine-readable ledgers with alignment flag, mapping quality, coordinate and CIGAR (Concise Idiosyncratic Gapped Alignment Report) representation.

Fifty-nine bounded local loci were assembled independently with MEGAHIT. All resulting contigs were aligned with `bwa-mem2 mem -a -Y` so that primary, secondary and supplementary contig placements were retained. External callers were used as algorithmic comparison: GRIDSS 2.13.2 used caller-compatible derived inputs with a retained filter-impact audit, and SvABA 1.2.0 used the corrected five-sample DLD1 invocation. The corrected SvABA run explicitly included all five DLD1 BAM files.

For every normalized external-caller breakend event, raw evidence was re-genotyped from the pre-existing, ungated full-window ledger for S01–S05. Caller fields, such as GRIDSS variant-fragment support, were reported as caller fields and were not re-labelled as independent-fragment counts. Re-genotyping and caller output derive from the same WGS data; they test mapping and algorithmic concordance, not independent biological replication.

## 7. Results

### 7.1 Coverage and full evidence inventory

All 18 samples completed the coverage and full-window inventory stages. Figure 1 shows canonical genome mean depth for every sample. Figure 2 shows the normalized chr11 trajectories for the five DLD1 comparison samples. The latter contains local reductions for S04/S05 after approximately chr11:114.604 Mb; this observation is not interpreted as a breakpoint or as quantitative confirmation of the qualitative fluorescence in situ hybridization (FISH) observation.

The raw evidence inventory retained 17,463–36,821 alignment records per sample in the discovery window. Mapping-quality-zero records comprised approximately 2.3–3.7% of the full-window record inventory, so excluding them a priori would have removed a non-trivial ambiguity track. Figure 3 documents this retention policy.

### 7.2 Assembly and placement audit

All 59 selected local assemblies completed. They produced 90 contigs, of which 32 had alternative genomic loci in the all-placement alignment audit. Figure 6 reports per-sample contigs with no mapped placement, one placement, or alternative loci. These results provide a transparent localization audit but no contig was promoted to a sequence-level breakpoint conclusion solely from its primary alignment.

Five unmapped mate read ends were recovered by complete source-BAM scans and incorporated as assembly input. Figure 5 records both the completed-locus count and this targeted rescue. The rescue preserves potentially informative sequence but does not increase the independent-fragment count by itself.

### 7.3 External caller and fixed-candidate evidence

The external-caller inventory contained 278 normalized GRIDSS events, four SvABA germline events and zero SvABA somatic events (Figure 7). Across the 282-event fixed universe, 94 events had two parsable endpoints; the remaining 188 were single-anchor or otherwise not evaluable as a two-endpoint geometry. These 188 records were retained in the audit rather than described as negative results.

The all-event GRIDSS DLD1 support matrix is shown in Figure 8. This panel displays an allocated caller support field, not a molecular count. Raw re-genotyping returned four case-only candidate patterns (Table 1 and Figure 9). The C2 pattern has three raw dual-endpoint QNAME fragments in S05, but none belongs to the standard-eligible subset; the other candidates have one to two `+dox` fragments. All controls have zero matching raw dual-endpoint fragments for these four patterns, but their standard local depths are only 14.7–23.1×. Consequently, the data do not support a 5% negative call in the controls.

| Candidate | Caller-defined endpoints (GRCh38) | Raw dual-endpoint fragments: S04 / S05 | Standard-eligible fragments: S04 / S05 | Control depth range (×) | Practical status |
| :--- | :---: | :---: | :---: | :---: | :---: |
| C1 | chr11:114564286; chr11:114564423 | 0 / 1 | 0 / 1 | 14.7–17.4 | Case-only raw candidate; insufficient confirmation |
| C2 | chr11:114607191; chr21:20648242 | 0 / 3 | 0 / 0 | 19.1–23.1 | Case-only raw candidate; evidence is not standard-eligible |
| C3 | chr11:114615021; chr11:114615035 | 1 / 1 | 0 / 1 | 19.3–22.1 | Case-only raw candidate; insufficient confirmation |
| C4 | chr11:114598286; chr3:69753603 | 0 / 1 | 0 / 1 | 16.5–17.4 | Case-only raw candidate; insufficient confirmation |

Figure 11 quantifies the fixed-candidate filter impact: 51,508 raw evidence rows were retained across S01–S05; 19,557 rows matched at least one caller anchor geometry and 31,951 did not. Rows not linked to a fixed caller geometry remain part of the full discovery evidence and were not deleted or interpreted as irrelevant.

## 8. Conclusions, Limitations and Next Steps

### 8.1 Detailed conclusions by decision question

| Decision question | Conclusion | Basis | Practical consequence |
| :--- | :---: | :---: | :--- |
| Is a Dox/clone-specific breakpoint confirmed? | No high-confidence breakpoint was established. | No candidate met the independent-fragment, control-sensitivity and sequence-level assembly gates together. | Do not use a WGS breakpoint as a definitive explanation of the FISH observation. |
| Are there focused hypotheses for follow-up? | Yes: C1–C4 are case-only raw-read candidates. | Exact two-endpoint raw evidence occurs only in S04/S05 under the fixed geometry rule. | Keep these four candidates as an explicit short list, not as confirmed structural variants. |
| Does coverage resolve the structure? | No. | S04/S05 coverage differs across portions of the chr11 window, but coverage has repeat and library-depth limitations. | Use coverage as context alongside, not instead of, junction evidence. |
| Is a unique contig placement demonstrated? | No. | 32/90 contigs have alternative loci. | Maintain alternative-placement annotation for any candidate carried forward. |

### 8.2 Practical reading with evidence sources

| Practical reading | What it adds beyond Section 8.1 | Evidence source |
| :--- | :---: | :--- |
| No confirmed Dox/clone-specific breakpoint | Lists the three separate gates that were not jointly satisfied. | `results/fixed_candidate_summary.tsv`; Figures 9, 10 and 12 |
| C1–C4 are focused candidate patterns | Gives exact GRCh38 coordinates and raw versus standard-eligible fragment counts. | `results/fixed_candidate_dld1_matrix.tsv`; Table 1; Figures 9 and 10 |
| Coverage is contextual evidence | Locates the sample-level chr11 trajectories and their normalization. | Figures 1 and 2; `figures/figure_manifest.json` |
| Alternative placement remains material | Provides the all-placement contig audit rather than a primary-only summary. | Figure 6; `results/contig_placement_summary.tsv` |

### 8.3 Limitations

This report does not quantify the FISH observation because a raw FISH image or tabular FISH result was not supplied in the analysis inputs. FISH and WGS are distinct evidence types; neither substitutes for the other.

The local control depth for C1–C4 is below the depth required for a three-fragment / 5% sensitivity context. Zero matching control fragments therefore means “not observed at available sensitivity,” not “biologically absent.” The four candidate patterns also lack a demonstrated sequence-level contig junction. Finally, all caller and re-genotyping evidence comes from the same WGS datasets and should not be described as independent biological confirmation.

### 8.4 Suggested confirmation

For a decision requiring a resolved structural event, prioritize the four listed candidates for an independent, breakpoint-spanning targeted assay. A targeted polymerase chain reaction (PCR) design should be based on the caller-defined two-endpoint geometry and assessed alongside the alternative-placement audit. If follow-up WGS is planned, higher local depth would improve the sensitivity context for comparisons against S01–S03.

## 9. Deliverable Files

| File or folder | Description |
| :--- | :--- |
| `P42_Vedita_18_WGS_1023.md` | This report; start here for the conclusions and file index. |
| `figures/` | Twelve numbered publication-style figures in Portable Network Graphics (PNG) and Portable Document Format (PDF), plus the figure manifest. |
| `results/fixed_candidate_summary.tsv` | One-row-per-candidate practical classification for all 282 fixed caller breakend events. |
| `results/fixed_candidate_dld1_matrix.tsv` | Five-sample raw and standard-eligible evidence matrix, local depth context and assembly-context fields. |
| `results/raw_evidence_match_audit.tsv.gz` | Machine-readable matched raw alignment-context audit, including QNAME, alignment class, mapping quality and matching rule. |
| `results/contig_placement_summary.tsv` | All-placement local-assembly contig summary, including alternative-locus information. |
| `methods/analysis_audit.json` | Fixed-candidate decision model, filter-impact counts and interpretation boundaries. |

Zhen Gao, PhD  
Principal Bioinformatics Scientist  
Athenomics
