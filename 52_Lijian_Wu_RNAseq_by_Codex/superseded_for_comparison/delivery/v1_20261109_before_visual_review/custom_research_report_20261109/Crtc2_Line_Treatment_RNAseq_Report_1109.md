# Bioinformatics Analysis Report: Control and Crtc2 Cell Lines under TNF, PGE2 and Adenosine

**Report Date:** 2026-11-09  
**Prepared by:** Zhen Gao, PhD, Principal Bioinformatics Scientist, Athenomics  
**Species:** Mus musculus (mouse; GRCm39, GENCODE M35)  
**Tissue / Cell:** Established mouse cell line (cell line name not specified by client)  
**Assay:** Bulk RNA sequencing (RNA-seq), 24 samples  
**Analysis Platform:** Linux HPC server  

## 1. Objectives

Compare the client-labelled control and Crtc2 cell lines under untreated (un), TNF, PGE2 and adenosine conditions, with three biological replicates per condition. CRTC2 denotes CREB-regulated transcription coactivator 2; CREB is cyclic adenosine monophosphate response element-binding protein. The modification in the Crtc2 line was not specified. The six requested comparisons use the second named group as baseline.

| Comparison | Numerator | Baseline |
| :--- | :---: | :---: |
| Crtc2-un vs control-un | Crtc2-un | control-un |
| Crtc2-TNF vs control-TNF | Crtc2-TNF | control-TNF |
| Crtc2-PGE2 vs control-PGE2 | Crtc2-PGE2 | control-PGE2 |
| Crtc2-adenosine vs control-adenosine | Crtc2-adenosine | control-adenosine |
| control-PGE2 vs control-un | control-PGE2 | control-un |
| control-adenosine vs control-un | control-adenosine | control-un |

## 2. Core Conclusions

This section explains the practical interpretation by experimental question. Section 3 summarises the measurements by topic. Section 8.2 adds the evidence files and inference boundaries for each conclusion; Section 8.3 describes limitations, and Section 8.4 identifies follow-up tests.

- **Quality and grouping.** All 24 libraries are retained. Unique alignment is 87.4–88.8%; the mean between-group / within-group distance in the first two principal components is 33.2. These checks support the expression analysis; they do not independently authenticate sample identity.

- **Crtc2 transcript.** Crtc2 messenger RNA (mRNA) is lower in the Crtc2 line in all four conditions (log2 fold changes -0.58, -0.93, -0.60, -0.59). The transcript remains detectable; protein loss and the nature of the modification are unresolved.

- **Untreated line difference.** The untreated lines differ in 1,151 genes. Lower inflammatory/interferon signatures and higher cell-cycle signatures are associated with the Crtc2 line; these are expression associations rather than a demonstrated mechanism.

- **TNF condition.** The largest line difference is under tumor necrosis factor (TNF): 5,906 genes. Inflammatory and interferon signatures are lower in the Crtc2 line under TNF. This comparison does not establish a difference in treatment responsiveness.

- **PGE2 and adenosine.** In control cells, prostaglandin E2 (PGE2) changes 2,175 genes and adenosine changes 254. Both show inflammatory/interferon enrichment and lower cell-cycle signatures; fewer adenosine genes cross the chosen gene-level threshold.

- **Practical next steps.** Confirm CRTC2 protein abundance and selected transcript changes in independent cultures. A formal line-by-treatment interaction analysis would address whether the two lines respond differently to treatment.

## 3. Key Findings

| Comparison | Up | Down | Total DEGs | DEGs with absolute log2FC ≥ 1 |
| :--- | :---: | :---: | :---: | :---: |
| Crtc2-un vs control-un | 461 | 690 | 1,151 | 246 |
| Crtc2-TNF vs control-TNF | 2,804 | 3,102 | 5,906 | 1,353 |
| Crtc2-PGE2 vs control-PGE2 | 799 | 1,307 | 2,106 | 539 |
| Crtc2-adenosine vs control-adenosine | 495 | 809 | 1,304 | 305 |
| control-PGE2 vs control-un | 1,278 | 897 | 2,175 | 625 |
| control-adenosine vs control-un | 170 | 84 | 254 | 13 |

A differentially expressed gene (DEG) meets adjusted p-value (padj) ≤ 0.05 and absolute log2 fold change (log2FC) ≥ 0.263, using shrunken fold changes. The Benjamini–Hochberg (BH) procedure controls the false discovery rate (FDR) within each comparison. The 0.263 cutoff is approximately 1.2-fold; it was explicitly selected for this analysis. The stricter twofold criterion from the form is also summarised above. Positive log2FC means higher expression in the numerator.

Quality control (QC) shows 28.5–65.1 million raw read pairs per sample, 96.0–97.1% total alignment and 0.03–0.06% ribosomal RNA (rRNA). Reverse library strandedness was inferred for all 24 samples.

Gene-level and pathway-level significance answer different questions: a treatment can have relatively few threshold-crossing genes while showing coordinated changes across gene sets. The complete pathway summary is provided in Section 7.3.

## 4. Sample Information

| Sample ID | Client sample name | Cell line | Treatment | Replicate | Raw read pairs | Estimated gene counts (rounded sum) | Detected genes |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| Control_1 | control-un-1 | control | un | 1 | 54,023,721 | 26,849,347 | 20,146 |
| Control_2 | control-un-2 | control | un | 2 | 50,135,520 | 26,113,022 | 20,167 |
| Control_3 | control-un-3 | control | un | 3 | 46,518,361 | 22,112,982 | 19,733 |
| Control_Ad_1 | control-adenosine-1 | control | adenosine | 1 | 38,221,188 | 20,238,228 | 19,245 |
| Control_Ad_2 | control-adenosine-2 | control | adenosine | 2 | 39,597,611 | 22,793,068 | 19,568 |
| Control_Ad_3 | control-adenosine-3 | control | adenosine | 3 | 32,738,797 | 23,984,817 | 18,635 |
| Control_PGE_1 | control-PEG2-1 | control | PGE2 | 1 | 33,541,400 | 20,359,577 | 19,671 |
| Control_PGE_2 | control-PEG2-2 | control | PGE2 | 2 | 45,530,162 | 22,842,077 | 19,743 |
| Control_PGE_3 | control-PEG2-3 | control | PGE2 | 3 | 50,076,161 | 27,142,559 | 20,201 |
| Control_TNFA_1 | control-TNF-1 | control | TNF | 1 | 43,181,663 | 20,239,906 | 19,771 |
| Control_TNFA_2 | control-TNF-2 | control | TNF | 2 | 28,532,014 | 17,056,583 | 19,348 |
| Control_TNFA_3 | control-TNF-3 | control | TNF | 3 | 33,245,777 | 22,102,489 | 20,302 |
| Crtc2_1 | Crtc2-un-1 | Crtc2 | un | 1 | 36,895,764 | 16,394,229 | 18,364 |
| Crtc2_2 | Crtc2-un-2 | Crtc2 | un | 2 | 47,989,941 | 25,023,342 | 19,524 |
| Crtc2_3 | Crtc2-un-3 | Crtc2 | un | 3 | 65,133,269 | 26,514,402 | 19,495 |
| Crtc2_Ad_1 | Crtc2-adenosine-1 | Crtc2 | adenosine | 1 | 52,955,882 | 26,315,086 | 19,563 |
| Crtc2_Ad_2 | Crtc2-adenosine-2 | Crtc2 | adenosine | 2 | 42,919,251 | 16,809,218 | 18,387 |
| Crtc2_Ad_3 | Crtc2-adenosine-3 | Crtc2 | adenosine | 3 | 34,909,566 | 17,782,254 | 18,008 |
| Crtc2_PGE_1 | Crtc2-PEG2-1 | Crtc2 | PGE2 | 1 | 46,250,803 | 20,890,211 | 19,200 |
| Crtc2_PGE_2 | Crtc2-PEG2-2 | Crtc2 | PGE2 | 2 | 46,140,077 | 17,911,902 | 18,782 |
| Crtc2_PGE_3 | Crtc2-PEG2-3 | Crtc2 | PGE2 | 3 | 49,575,310 | 21,579,191 | 19,215 |
| Crtc2_TNFA_1 | Crtc2-TNF-1 | Crtc2 | TNF | 1 | 50,290,839 | 18,277,584 | 19,402 |
| Crtc2_TNFA_2 | Crtc2-TNF-2 | Crtc2 | TNF | 2 | 48,985,421 | 20,532,208 | 19,804 |
| Crtc2_TNFA_3 | Crtc2-TNF-3 | Crtc2 | TNF | 3 | 43,078,618 | 23,141,234 | 20,017 |

Client sample names retain the form spelling PEG2; PGE2 is used for analysis labels, as written in the requested comparisons. FASTQ-to-client sample mapping was confirmed before analysis. The form specifies established mouse cells, snap-frozen. The specific cell line and Crtc2 modification were not provided. The assigned-count sum is based on rounded Salmon gene estimates and is not a direct raw-fragment count.

## 5. Analysis Rationale and Decision Criteria

| Step | Criterion | Reason and scope |
| :--- | :---: | :---: |
| Reference | GRCm39 primary assembly; GENCODE M35 | Matches the stated species and the reference used for quantification. |
| Gene-name filter | 57,132 → 26,584 genes | Exclude ribosomal-protein names and specified mitochondrial/noncoding/predicted-name patterns. This is a name-based restriction, not a complete biotype classification; excluded loci remain in the supplied unfiltered matrices. |
| Low-expression filter | 26,584 → 12,327 genes in the model union | For each comparison, ≥10 rounded estimated counts in ≥3 of the 6 compared samples OR ≥8 of all 24 samples. This retains condition-restricted expression without requiring expression in most groups. |
| Model | One DESeq2 ~ Group model, eight levels | Six predefined contrasts share the fitted model. Biological samples, not reads, are the replicate units; no samples were removed. |
| Multiplicity | independentFiltering=FALSE; BH within the eligible set | Each comparison uses its own explicitly defined testable genes. |
| Effect size | ashr shrinkage; absolute log2FC ≥ 0.263 | Stabilise effect estimates and apply the confirmed effect-size threshold alongside padj. |
| Over-representation analysis (ORA) | UP / DOWN / ALL; comparison-specific mapped background | Compare the selected gene set with genes that were testable and could be mapped to identifiers. |
| Gene set enrichment analysis (GSEA) | Unshrunken DESeq2 Wald statistic; seed 42 | Assess coordinated ranked expression shifts without restricting input to DEGs. Enrichment is an association, not an activity assay. |

## 6. Methods

| Step | Software / reference | Parameters |
| :--- | :---: | :---: |
| Pipeline | nf-core/rnaseq 3.15.1 | STAR alignment with Salmon quantification; paired-end reads; strandedness inferred automatically |
| Alignment | STAR; GRCm39 / GENCODE M35 | --twopassMode None --outFilterMultimapNmax 3 |
| Preprocessing and QC | Trim Galore!, FastQC, RSeQC, Qualimap, Picard, MultiQC | Pipeline-generated trimming, alignment and library QC summaries |
| Differential expression | DESeq2 1.50.2; ashr 2.2.63 | Wald test; ~ Group; explicit per-comparison eligibility and BH correction; shrunken log2FC |
| ORA / GSEA | clusterProfiler 4.18.4; fgsea 1.36.2 | Gene Ontology (GO) biological process (BP), molecular function (MF), cellular component (CC); Kyoto Encyclopedia of Genes and Genomes (KEGG); Molecular Signatures Database (MSigDB) Hallmark |
| Mouse identifier mapping | org.Mm.eg.db 3.22.0 | Entrez identifier mapping; GSEA duplicate identifiers represented by the largest absolute Wald statistic |

The supplied count matrix contains unrounded Salmon gene-level estimated counts. The expression filter and DESeq2 use integer-rounded estimates. Transcripts per million (TPM) are supplied separately. Gene-type and low-expression filters apply to statistical testing, not to the unfiltered delivered count/TPM matrices. ORA and GSEA tables contain terms passing their output thresholds; missing terms must not be read as demonstrated absence of a biological process.

## 7. Results

### 7.1 Principal Component Analysis (PCA)

The first principal component (PC1) explains 93% of the variance and PC2 explains 4%. PCA uses variance-stabilised expression of the top 500 variable model genes. The overview includes 24 points and 24 sample labels. Mean within-group distance is 0.62 and mean between-group distance is 20.60 in PC1–PC2 space (ratio 33.2). Coordinates are in `DE_PCA_Results/PCA_coordinates_overview.csv`; the figure is `DE_PCA_Results/PCA_AllSamples_Overview.pdf`.

| Group | Mean PC1 | Mean PC2 |
| :--- | :---: | :---: |
| Crtc2-PGE2 | -9.53 | 0.31 |
| Crtc2-TNF | 19.06 | 6.32 |
| Crtc2-adenosine | -13.51 | 2.61 |
| Crtc2-un | -13.72 | 2.34 |
| control-PGE2 | -0.48 | -5.29 |
| control-TNF | 36.86 | -2.42 |
| control-adenosine | -8.11 | -1.81 |
| control-un | -10.56 | -2.06 |

TNF groups occupy the high-PC1 region, making treatment the dominant separation in this view. PC2 separates the two cell lines. These patterns support broad biological structure but do not substitute for a sample-identity assay.

Nearest-neighbour checks found 1 sample outside its own group: Crtc2_1 → Crtc2_Ad_3. This involves the closely overlapping Crtc2 untreated and adenosine groups; the global two-component view does not resolve these conditions well. All other samples have a nearest neighbour in their own group. This isolated overlap is insufficient evidence of mislabelling, and all samples were retained.

### 7.2 Differential Expression

| Comparison | Genes tested | DEGs | DEGs / tested | Crtc2 log2FC | Crtc2 padj |
| :--- | :---: | :---: | :---: | :---: | :---: |
| Crtc2-un vs control-un | 11,883 | 1,151 | 9.7% | -0.58 | 1.24e-25 |
| Crtc2-TNF vs control-TNF | 12,172 | 5,906 | 48.5% | -0.93 | 4.81e-83 |
| Crtc2-PGE2 vs control-PGE2 | 11,878 | 2,106 | 17.7% | -0.60 | 1.11e-29 |
| Crtc2-adenosine vs control-adenosine | 11,823 | 1,304 | 11.0% | -0.59 | 3.18e-26 |
| control-PGE2 vs control-un | 11,937 | 2,175 | 18.2% | +0.38 | 9.26e-14 |
| control-adenosine vs control-un | 11,872 | 254 | 2.1% | +0.01 | 0.897 |

The TNF line comparison has the largest DEG count. The two control-treatment comparisons show a much larger threshold-crossing response to PGE2 than to adenosine. These counts quantify the six specified comparisons; comparing DEG totals does not itself test an interaction.

| Condition | Control mean Crtc2 normalised count | Crtc2-line mean | Crtc2 line / control |
| :--- | :---: | :---: | :---: |
| un | 1376.6 | 919.9 | 0.67 |
| TNF | 4078.1 | 2141.4 | 0.53 |
| PGE2 | 1816.7 | 1195.1 | 0.66 |
| adenosine | 1433.1 | 942.8 | 0.66 |

These target-gene group means use library-size normalisation estimated from the full count matrix, as provided in `Reads/Crtc2_expression_per_sample.csv`. The DE log2FC values above use the filtered DESeq2 model and shrinkage, so the two summaries need not be algebraically identical.

Every comparison supplies a full tested-gene table, significant subset, volcano plot, and top-50 DEG heatmap. Genes outside a comparison’s eligible set were not tested there, rather than demonstrated unchanged.

### 7.3 Pathway Enrichment

Hallmark normalised enrichment score (NES): positive means enrichment toward the numerator; negative toward the baseline. “n.s.” means the gene set did not appear among the significant outputs (padj < 0.05), not absence of the process.

| Comparison | TNFA SIGNALING VIA NFKB | INFLAMMATORY RESPONSE | INTERFERON GAMMA RESPONSE | E2F TARGETS | G2M CHECKPOINT |
| :--- | :---: | :---: | :---: | :---: | :---: |
| Crtc2-un vs control-un | -1.57 | -2.15 | -1.65 | +2.00 | +1.86 |
| Crtc2-TNF vs control-TNF | -2.10 | -1.62 | -1.96 | n.s. | n.s. |
| Crtc2-PGE2 vs control-PGE2 | -2.21 | -2.22 | -2.12 | +2.10 | +1.97 |
| Crtc2-adenosine vs control-adenosine | -1.81 | -2.24 | -1.98 | +2.36 | +2.39 |
| control-PGE2 vs control-un | +2.59 | +2.19 | +2.34 | -3.31 | -3.15 |
| control-adenosine vs control-un | +2.75 | +2.17 | +2.46 | -3.79 | -3.53 |

Inflammatory and interferon-related signatures are lower in the Crtc2 line in the four matched-condition comparisons. Higher cell-cycle signatures are seen in the untreated, PGE2 and adenosine line comparisons; the TNF condition does not show the same significant E2F/G2M pattern. In control cells, both PGE2 and adenosine show the opposite direction for these themes. These results motivate functional follow-up but do not prove direct CRTC2 regulation, pathway activity, or cell proliferation.

| Comparison | GO_BP_UP | GO_BP_DOWN | KEGG_UP | KEGG_DOWN | GSEA_GO_BP | GSEA_KEGG | GSEA_Hallmark |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| Crtc2-un vs control-un | 113 | 538 | 8 | 17 | 295 | 12 | 16 |
| Crtc2-TNF vs control-TNF | 137 | 516 | 14 | 17 | 73 | 23 | 14 |
| Crtc2-PGE2 vs control-PGE2 | 177 | 1067 | 12 | 38 | 476 | 32 | 19 |
| Crtc2-adenosine vs control-adenosine | 186 | 604 | 27 | 28 | 473 | 4 | 18 |
| control-PGE2 vs control-un | 1024 | 253 | 43 | 12 | 1172 | 107 | 28 |
| control-adenosine vs control-un | 273 | 203 | 20 | 4 | 1033 | 101 | 25 |

Counts are significant terms at adjusted p-value < 0.05; overlapping terms are not independent biological discoveries. The Hallmark display above is a selected thematic view; all returned terms and the full 90-analysis summary are provided in `Enrichment/`. ORA UP / DOWN / ALL are separate tests, with the mapped testable genes as background.

## 8. Conclusions, Limitations and Next Steps

### 8.1 Per-comparison interpretation

| Comparison | Meaning | Boundary |
| :--- | :---: | :---: |
| Crtc2-un vs control-un | Baseline differences include inflammatory and cell-cycle expression themes. | Expression association within this comparison; no formal interaction test. |
| Crtc2-TNF vs control-TNF | Largest matched-condition line difference; lower inflammatory and interferon signatures in Crtc2. | Expression association within this comparison; no formal interaction test. |
| Crtc2-PGE2 vs control-PGE2 | Line difference under PGE2; it does not alone establish a difference in PGE2 responsiveness. | Expression association within this comparison; no formal interaction test. |
| Crtc2-adenosine vs control-adenosine | Line difference under adenosine; the Crtc2 untreated and adenosine groups overlap in PC1–PC2. | Expression association within this comparison; no formal interaction test. |
| control-PGE2 vs control-un | Broad control-cell PGE2-associated expression change, including inflammatory and cell-cycle themes. | Expression association within this comparison; no formal interaction test. |
| control-adenosine vs control-un | Fewer gene-level DEGs than PGE2, with coordinated pathway shifts despite small gene-level effects. | Expression association within this comparison; no formal interaction test. |

### 8.2 Practical reading

The following evidence map gives the data location and inference boundary for each conclusion in Section 2.

| Question | Measured evidence / delivered source | Interpretation and assumptions |
| :--- | :---: | :---: |
| Quality and grouping | Section 7.1; QC/multiqc/; DE_PCA_Results/PCA_coordinates_overview.csv | None; no identity assay was performed. |
| Crtc2 transcript | Section 7.2; DE_PCA_Results/<comparison>/<comparison>_ALL_genes.csv; row Crtc2 | Protein and functional effects require independent confirmation. |
| Untreated line difference | Section 7.2 and 7.3; DE_PCA_Results/Crtc2-un_vs_control-un/; Enrichment/Crtc2-un_vs_control-un/GSEA/ | An association with line identity; not proof of direct CRTC2 regulation. |
| TNF condition | Section 7.2 and 7.3; DE_PCA_Results/Crtc2-TNF_vs_control-TNF/; Enrichment/Crtc2-TNF_vs_control-TNF/GSEA/ | A line-by-treatment interaction was not tested. |
| PGE2 and adenosine | Section 7.2 and 7.3; DE_PCA_Results/control-PGE2_vs_control-un/; DE_PCA_Results/control-adenosine_vs_control-un/; Enrichment/ | Pathway enrichment does not measure pathway activity or establish equivalent effect sizes. |
| Practical next steps | Section 8.4; Section 8.4 | These are proposed follow-up tests, not completed validation. |

### 8.3 Limitations

- Three biological replicates per condition limit power for small or variable changes. A nonsignificant or untested gene is not shown to be unchanged.
- The specific cell line and the nature of the Crtc2 modification were not supplied. Reduced Crtc2 mRNA does not demonstrate loss of CRTC2 protein or function.
- Only six requested contrasts were tested. There is no formal line-by-treatment interaction test or delivered within-line TNF-response contrast. PCA overlap and differing DEG totals do not establish equivalent responses or loss of responsiveness.
- Name-based filtering restricts the statistical gene universe and is not a complete protein-coding biotype annotation. Low-expression and identifier-mapping filters also limit coverage. Filtered genes remain available in the full expression matrices.
- Pathway terms overlap and reflect gene-expression associations. Enrichment does not prove pathway activation, altered proliferation, or a direct molecular mechanism.
- Alignment and expression QC do not independently authenticate sample identity. The overlap of Crtc2 untreated and adenosine samples is reported explicitly rather than treated as proof of either mislabelling or equivalence.

### 8.4 Suggested confirmation

| Open question | Follow-up |
| :--- | :---: |
| Is CRTC2 protein reduced? | Western blot comparing control and Crtc2 cells under matched conditions. |
| Do the transcript changes replicate? | Quantitative polymerase chain reaction (qPCR) for Crtc2 and selected DEGs in independent cultures. |
| Does treatment responsiveness differ between lines? | A prespecified line-by-treatment interaction analysis on the existing samples; this would be an additional statistical analysis, not independent biological validation. |
| Do pathway themes have functional consequences? | Independent proliferation measurements and cytokine-protein assays selected for the biological question. |

## 9. Deliverable Files

| Path | Content |
| :--- | :---: |
| `Reads/All_sample_gene_counts.tsv` | Unfiltered Salmon gene-level estimated counts, 24 samples |
| `Reads/All_sample_gene_tpm.tsv` | Unfiltered gene-level TPM, 24 samples |
| `Reads/library_metrics_per_sample.csv` | Raw read pairs, estimated-count sums and detected genes |
| `Reads/Crtc2_expression_per_sample.csv` | Target-gene normalised counts for 24 samples |
| `QC/multiqc/` | Pipeline MultiQC HTML report and supporting data |
| `DE_PCA_Results/PCA_AllSamples_Overview.pdf` | 24-sample PCA |
| `DE_PCA_Results/PCA_coordinates_overview.csv` | 24 sample coordinates and grouping |
| `DE_PCA_Results/SampleDistance_Heatmap.pdf` | Sample distances with group annotation |
| `DE_PCA_Results/normalized_counts_all_samples.csv` | Filtered-model normalised counts |
| `DE_PCA_Results/DEG_summary_all_comparisons.csv` | Six-comparison gene counts at both effect-size thresholds |
| `DE_PCA_Results/<comparison>/` | Full and significant tables, volcano and top-50 heatmap for each of six comparisons |
| `Enrichment/<comparison>/{GO,KEGG,GSEA}/` | Directional ORA and ranked GSEA tables and PDF plots |
| `Enrichment/Enrichment_summary_all_comparisons.csv` | 90-analysis status and significant-term counts |
| `Gene_filtering/` | Full-gene and six per-comparison eligibility audits |
| `Claim_evidence_map.csv` | Measured facts, interpretation boundaries and source files for Section 2 |
| `mouse_Gene_annotation_20260628.xlsx` | Gene annotation reference |
| `Crtc2_Line_Treatment_RNAseq_Report_1109.md` | Report and reading entry point |

---

*Zhen Gao, PhD, Principal Bioinformatics Scientist, Athenomics*
