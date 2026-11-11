# Bioinformatics Analysis Report: Control and Crtc2 Cell Lines under TNF, PGE2 and Adenosine

**Report Date:** 2026-11-08  
**Prepared by:** Zhen Gao, PhD, Principal Bioinformatics Scientist, Athenomics  
**Species:** Mus musculus (mouse; reference genome GRCm39, GENCODE M35 annotation)  
**Tissue / Cell:** Established mouse cell line (cell line name not specified by client)  
**Assay:** Bulk RNA sequencing (RNA-seq), poly-A-selected mRNA, 24 samples  
**Analysis Platform:** Linux HPC server  

---

## 1. Objectives

The study compares two mouse cell lines, labelled **control** and **Crtc2** by the client, each left untreated (**un**) or treated with **TNF** (tumor necrosis factor), **PGE2** (prostaglandin E2) or **adenosine**; three biological replicates per condition, 24 samples in total. Crtc2 encodes CRTC2 (CREB-regulated transcription coactivator 2), a co-activator of the transcription factor CREB (cAMP response element-binding protein; cAMP = cyclic adenosine monophosphate). The submission form does not state how the Crtc2 line was modified, so this report refers to it as the **Crtc2 line** and reports the Crtc2 mRNA level that was measured.

The six comparisons requested on the submission form (the second group named is the baseline):

| # | Comparison | Question it answers |
| :--- | :---: | :---: |
| 1 | Crtc2-un vs control-un | How the Crtc2 line differs from control without treatment |
| 2 | Crtc2-TNF vs control-TNF | How the Crtc2 line differs from control under TNF |
| 3 | Crtc2-PGE2 vs control-PGE2 | How the Crtc2 line differs from control under PGE2 |
| 4 | Crtc2-adenosine vs control-adenosine | How the Crtc2 line differs from control under adenosine |
| 5 | control-PGE2 vs control-un | The PGE2 response of control cells |
| 6 | control-adenosine vs control-un | The adenosine response of control cells |

## 2. Core Conclusions

This section gives the answer for each part of the experiment in plain language. Section 3 lists what was measured topic by topic; Section 8.2 gives **more detail on each conclusion below and names the section and delivered file that supports it**; Section 8.3 states what these data cannot tell you; Section 8.4 lists the experiment that would settle each open question.

- **The data are of good quality and the sample labels are consistent.** All 24 libraries passed quality control (87.4–88.8% of reads uniquely mapped), all were sequenced together on one flowcell lane, and replicates of the same condition are 33.2 times closer to each other than to other conditions. No sample needs to be excluded.

- **The Crtc2 line has 33–47% less Crtc2 mRNA than control, in every condition** (Crtc2 line / control = 0.53–0.67). The transcript is reduced, not absent. RNA-seq measures mRNA only, so whether CRTC2 protein is lost should be confirmed by Western blot (Section 8.4).

- **Without treatment, the Crtc2 line already differs from control in 1,151 genes, and the difference runs opposite to the PGE2 response.** Classical CREB target genes (Nr4a2, Crem, Pde4b, Pde4d, Sik1) are lower in the Crtc2 line; of the 740 genes changed both by the Crtc2 line and by PGE2, 91% change in opposite directions.

- **PGE2 is a strong stimulus in control cells (2,174 genes change), and the Crtc2 line reproduces only about 58% of that response.** In control cells PGE2 raises CREB targets and inflammatory genes (Pde4d +2.2, Il6 +2.7 log2 fold change) and lowers cell-cycle genes. Under PGE2 the Crtc2 line differs from control in 2,105 genes; 66% of them are genes that PGE2 changes in control cells, and 97% of those move in the opposite direction.

- **Adenosine has a mild effect on control cells (254 genes, only 13 of them changing 2-fold or more) and almost none on the Crtc2 line.** In the Crtc2 line, adenosine-treated and untreated samples are practically indistinguishable (Section 7.1); on the control-cell adenosine genes the Crtc2 line reproduces about 32% of the response.

- **TNF is the strongest stimulus in the experiment, and the largest Crtc2-line difference is seen under TNF (5,905 genes).** The Crtc2 line still responds to TNF, but less (about 74% of the control response), and TNF-induced inflammatory genes such as Il6 (-3.4 log2 fold change) and Cxcl10 are clearly lower than in control.

- **What we suggest next.** Confirm CRTC2 protein loss by Western blot, and confirm a small panel of the key genes (Pde4d, Nr4a2, Il6) by qPCR (quantitative PCR) in independent cultures.

## 3. Key Findings

Section 2 gave the conclusion for each part of the experiment. This section reports the measured results by topic.

**Differentially expressed genes (DEGs) per comparison**

| Comparison | Up-regulated | Down-regulated | Total DEGs | DEGs with at least 2-fold change |
| :--- | :---: | :---: | :---: | :---: |
| Crtc2-un vs control-un | 461 | 690 | **1,151** | 246 (56 up / 190 down) |
| Crtc2-TNF vs control-TNF | 2,804 | 3,101 | **5,905** | 1,353 (556 up / 797 down) |
| Crtc2-PGE2 vs control-PGE2 | 798 | 1,307 | **2,105** | 539 (102 up / 437 down) |
| Crtc2-adenosine vs control-adenosine | 495 | 809 | **1,304** | 305 (53 up / 252 down) |
| control-PGE2 vs control-un | 1,278 | 896 | **2,174** | 625 (436 up / 189 down) |
| control-adenosine vs control-un | 170 | 84 | **254** | 13 (11 up / 2 down) |

A DEG is a gene with adjusted p-value (padj; Benjamini-Hochberg false discovery rate, FDR) of 0.05 or less and a change of at least 1.2-fold (absolute log2 fold change, log2FC, of at least 0.263). The last column applies the stricter 2-fold cut-off (absolute log2FC of at least 1) shown as the default on the submission form. The 1.2-fold cut-off is used as the primary criterion because many consistent, biologically relevant changes in this data set are smaller than 2-fold; Crtc2 itself (log2FC about -0.6) would not pass a 2-fold cut-off. A positive log2FC means higher expression in the first group named.

**Crtc2 mRNA.** Lower in the Crtc2 line in all four conditions (log2FC -0.58, -0.93, -0.60, -0.59 for comparisons 1–4; all significant).

**How much of each treatment response the Crtc2 line keeps** (descriptive; Section 7.4):

| Treatment | Genes used | Crtc2-line response as a fraction of the control response |
| :--- | :---: | :---: |
| PGE2 | 2,174 | 58% |
| adenosine | 254 | 32% |
| TNF | 3,950 | 74% |

**Pathway themes.** In control cells, PGE2 and adenosine both raise inflammatory and interferon gene programmes and lower cell-cycle programmes. In the Crtc2 line the same programmes move the opposite way: inflammatory and interferon genes are lower and cell-cycle genes higher than in control, at baseline and under every treatment (Section 7.3).

**Data quality.** 28.5–65.1 million read pairs per sample (paired-end, 150 bp), 96.0–97.1% of reads aligned, 87.4–88.8% uniquely aligned, rRNA (ribosomal RNA) 0.03–0.06%. Library strandedness was inferred from the data as reverse for all 24 samples.

## 4. Sample Information

| Sample ID (this report) | Client sample name | Cell line | Treatment | Replicate | Read pairs | Reads assigned to genes | Genes detected |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| control-un-1 | control-un-1 | control | un | 1 | 54,023,721 | 26,849,344 | 20,143 |
| control-un-2 | control-un-2 | control | un | 2 | 50,135,520 | 26,113,026 | 20,169 |
| control-un-3 | control-un-3 | control | un | 3 | 46,518,361 | 22,112,983 | 19,733 |
| control-PGE2-1 | control-PEG2-1 | control | PGE2 | 1 | 33,541,400 | 20,359,580 | 19,668 |
| control-PGE2-2 | control-PEG2-2 | control | PGE2 | 2 | 45,530,162 | 22,842,082 | 19,745 |
| control-PGE2-3 | control-PEG2-3 | control | PGE2 | 3 | 50,076,161 | 27,142,570 | 20,200 |
| control-TNF-1 | control-TNF-1 | control | TNF | 1 | 43,181,663 | 20,239,913 | 19,771 |
| control-TNF-2 | control-TNF-2 | control | TNF | 2 | 28,532,014 | 17,056,576 | 19,350 |
| control-TNF-3 | control-TNF-3 | control | TNF | 3 | 33,245,777 | 22,102,491 | 20,302 |
| control-adenosine-1 | control-adenosine-1 | control | adenosine | 1 | 38,221,188 | 20,238,224 | 19,244 |
| control-adenosine-2 | control-adenosine-2 | control | adenosine | 2 | 39,597,611 | 22,793,061 | 19,570 |
| control-adenosine-3 | control-adenosine-3 | control | adenosine | 3 | 32,738,797 | 23,984,802 | 18,633 |
| Crtc2-un-1 | Crtc2-un-1 | Crtc2 | un | 1 | 36,895,764 | 16,394,217 | 18,367 |
| Crtc2-un-2 | Crtc2-un-2 | Crtc2 | un | 2 | 47,989,941 | 25,023,348 | 19,525 |
| Crtc2-un-3 | Crtc2-un-3 | Crtc2 | un | 3 | 65,133,269 | 26,514,399 | 19,494 |
| Crtc2-PGE2-1 | Crtc2-PEG2-1 | Crtc2 | PGE2 | 1 | 46,250,803 | 20,890,211 | 19,196 |
| Crtc2-PGE2-2 | Crtc2-PEG2-2 | Crtc2 | PGE2 | 2 | 46,140,077 | 17,911,908 | 18,783 |
| Crtc2-PGE2-3 | Crtc2-PEG2-3 | Crtc2 | PGE2 | 3 | 49,575,310 | 21,579,181 | 19,215 |
| Crtc2-TNF-1 | Crtc2-TNF-1 | Crtc2 | TNF | 1 | 50,290,839 | 18,277,585 | 19,402 |
| Crtc2-TNF-2 | Crtc2-TNF-2 | Crtc2 | TNF | 2 | 48,985,421 | 20,532,214 | 19,807 |
| Crtc2-TNF-3 | Crtc2-TNF-3 | Crtc2 | TNF | 3 | 43,078,618 | 23,141,222 | 20,013 |
| Crtc2-adenosine-1 | Crtc2-adenosine-1 | Crtc2 | adenosine | 1 | 52,955,882 | 26,315,087 | 19,563 |
| Crtc2-adenosine-2 | Crtc2-adenosine-2 | Crtc2 | adenosine | 2 | 42,919,251 | 16,809,205 | 18,384 |
| Crtc2-adenosine-3 | Crtc2-adenosine-3 | Crtc2 | adenosine | 3 | 34,909,566 | 17,782,253 | 18,009 |

The client sample names list the PGE2 samples as "PEG2" (for example control-PEG2-1), while the comparison list on the same form writes "pge2"; this report uses **PGE2** throughout and keeps the client name alongside. All samples are from an established mouse cell line ("Normal only", snap-frozen, per the form). All 24 libraries were sequenced on the same flowcell lane (2533TCLT4_3), so sequencing run cannot be confounded with cell line or treatment.

## 5. Analysis Rationale and Decision Criteria

| Step | Choice | Why |
| :--- | :---: | :---: |
| Reference | GRCm39 primary assembly, GENCODE M35 | Current mouse reference genome and annotation |
| Strandedness | Inferred from the data | Avoids assuming the library chemistry; a wrong setting would distort gene counts |
| Statistical model | One DESeq2 model, `~ Group`, all 24 samples (8 groups) | Single sequencing lane and balanced design; all six comparisons share one variance estimate from 24 samples, which is more precise than fitting each pair separately |
| Gene-type filter | 57,132 -> 26,584 genes | Removes ribosomal-protein, mitochondrial and non-coding/predicted loci, which add noise and dilute multiple-testing correction |
| Low-expression filter | 26,584 -> 12,327 genes (11,823–12,172 tested per comparison) | A gene is tested in a comparison if it has at least 10 counts (Salmon estimated counts rounded to integers, as used by DESeq2) in at least 3 of the 6 samples being compared, or in at least 30% of all samples. This keeps genes that are switched on by one treatment only, which a rule requiring expression across most of the 24 samples would discard |
| Significance | padj of 0.05 or less and at least 1.2-fold change | False-discovery control plus a minimum effect size; counts at the 2-fold cut-off are also reported |
| Fold-change shrinkage | `ashr` | Stabilises fold changes of low-count genes |
| Over-representation analysis | Up- and down-regulated genes tested separately | Testing both directions together loses power when many genes respond |
| GSEA ranking | DESeq2 Wald statistic | Uses every tested gene without a cut-off; the Wald statistic weights genes evenly enough that a few extreme genes cannot dominate the result |

## 6. Methods

| Step | Tool | Key parameters |
| :--- | :---: | :---: |
| Pipeline | nf-core/rnaseq 3.15.1 | `--aligner star_salmon --gencode` |
| Trimming | Trim Galore! | nf-core defaults |
| Alignment | STAR | GRCm39 index; `--outFilterMultimapNmax 3` |
| Quantification | Salmon | Gene-level estimated counts and TPM (transcripts per million) |
| Quality control | FastQC, RSeQC, Qualimap, Picard, MultiQC | Strandedness inferred by RSeQC and Salmon |
| Differential expression | DESeq2 | `~ Group`; Wald test; per-comparison gene set as in Section 5 with Benjamini-Hochberg correction within it; `lfcShrink(type = "ashr")` |
| Over-representation | clusterProfiler | GO (Gene Ontology) biological process (BP), molecular function (MF) and cellular component (CC), and KEGG (Kyoto Encyclopedia of Genes and Genomes); background = genes tested in that comparison |
| GSEA (gene set enrichment analysis) | clusterProfiler / fgsea | Ranked by the DESeq2 Wald statistic; GO BP, KEGG, MSigDB (Molecular Signatures Database) Hallmark gene sets (mouse orthologs) |

The count matrix contains Salmon **gene-level estimated counts**, not raw read counts (`Reads/All_sample_gene_counts.tsv`, delivered unrounded); DESeq2 and the low-expression filter use them rounded to integers. `Reads/All_sample_gene_tpm.tsv` contains TPM values.

## 7. Results

### 7.1 Principal Component Analysis (PCA)

The first principal component (PC1) explains 93% of the variance and the second (PC2) 4% (`DE_PCA_Results/PCA_AllSamples_Overview.pdf`; coordinates in `DE_PCA_Results/PCA_coordinates_overview.csv`). Replicates of the same condition lie on average 0.62 units apart and samples of different conditions 20.58 units apart, a ratio of 33.2: the conditions are far more distinct than the replicate-to-replicate variation.

**PC1 is the TNF axis.** The TNF-treated samples sit at one end (mean PC1 27.9) and all other samples at the other (untreated -12.1, PGE2 -5.0, adenosine -10.8). Along this axis the Crtc2-TNF group lies between control-TNF and the untreated groups: control-TNF is 47.4 units from control-un, Crtc2-TNF 33.0 units from Crtc2-un, so the TNF shift is smaller in the Crtc2 line.

Among the non-TNF groups the same pattern holds for PGE2 and adenosine. PGE2 moves control cells 10.6 units from control-un but the Crtc2 line only 4.6 units from Crtc2-un; adenosine moves control cells 2.4 units but the Crtc2 line only 0.34 units, which is no further apart than replicates of one condition. The Crtc2-PGE2 group sits only 2.6 units from untreated control cells.

**PC2 is the cell-line axis.** Every Crtc2-line sample has a higher PC2 value (0.1 to 6.4) than every control sample (-5.5 to -1.7), in all four conditions: the second-largest source of variation in the data set is the difference between the two lines.

Each sample's nearest neighbour in the PCA was checked. 1 of 24 samples sits closest to a sample of another condition (Crtc2-un-1 -> Crtc2-adenosine-3). This involves only the Crtc2-un and Crtc2-adenosine groups, which are indistinguishable because adenosine has no measurable effect on the Crtc2 line (above); it is not a sign of a sample mix-up. Every other sample sits closest to a replicate of its own condition.

| Pair of conditions | PCA centroid distance |
| :--- | :---: |
| control-un vs control-TNF | 47.36 |
| Crtc2-un vs Crtc2-TNF | 32.96 |
| control-un vs control-PGE2 | 10.55 |
| Crtc2-un vs Crtc2-PGE2 | 4.64 |
| control-un vs control-adenosine | 2.45 |
| Crtc2-un vs Crtc2-adenosine | 0.34 |
| control-un vs Crtc2-un | 5.42 |
| control-TNF vs Crtc2-TNF | 19.82 |
| control-PGE2 vs Crtc2-PGE2 | 10.63 |
| control-adenosine vs Crtc2-adenosine | 6.98 |
| control-un vs Crtc2-PGE2 | 2.58 |

### 7.2 Differential Expression

| Comparison | Genes tested | Up | Down | Total DEGs | % of tested genes |
| :--- | :---: | :---: | :---: | :---: | :---: |
| Crtc2-un vs control-un | 11,883 | 461 | 690 | **1,151** | 9.7% |
| Crtc2-TNF vs control-TNF | 12,172 | 2,804 | 3,101 | **5,905** | 48.5% |
| Crtc2-PGE2 vs control-PGE2 | 11,878 | 798 | 1,307 | **2,105** | 17.7% |
| Crtc2-adenosine vs control-adenosine | 11,823 | 495 | 809 | **1,304** | 11.0% |
| control-PGE2 vs control-un | 11,937 | 1,278 | 896 | **2,174** | 18.2% |
| control-adenosine vs control-un | 11,872 | 170 | 84 | **254** | 2.1% |

Each comparison folder in `DE_PCA_Results/` holds the full table of tested genes (fold change, Wald statistic, p-value, padj and normalised counts of the six samples), the significant subset, a volcano plot and a heatmap of the top 50 DEGs.

**Crtc2 mRNA** (DESeq2-normalised counts, mean of three replicates):

| Treatment | control | Crtc2 line | Crtc2 line / control |
| :--- | :---: | :---: | :---: |
| untreated | 1,377 | 920 | 0.67 |
| TNF | 4,078 | 2,141 | 0.53 |
| PGE2 | 1,817 | 1,195 | 0.66 |
| adenosine | 1,433 | 943 | 0.66 |

Crtc2 mRNA is higher with TNF in both lines (group means 3.0-fold in control, 2.3-fold in the Crtc2 line; TNF vs untreated was not a tested comparison). In every condition the Crtc2 line carries 33–47% less Crtc2 mRNA than control.

### 7.3 Pathway Enrichment

**Hallmark gene sets (GSEA normalised enrichment score, NES; positive = higher in the first group named; n.s. = not significant)**

| Comparison | TNFA SIGNALING VIA NFKB | INFLAMMATORY RESPONSE | INTERFERON GAMMA RESPONSE | E2F TARGETS | G2M CHECKPOINT |
| :--- | :---: | :---: | :---: | :---: | :---: |
| Crtc2-un vs control-un | -1.56 | -2.15 | -1.65 | +2.00 | +1.86 |
| Crtc2-TNF vs control-TNF | -2.10 | -1.62 | -1.96 | n.s. | n.s. |
| Crtc2-PGE2 vs control-PGE2 | -2.21 | -2.22 | -2.12 | +2.10 | +1.97 |
| Crtc2-adenosine vs control-adenosine | -1.81 | -2.24 | -1.98 | +2.36 | +2.39 |
| control-PGE2 vs control-un | +2.59 | +2.19 | +2.35 | -3.31 | -3.15 |
| control-adenosine vs control-un | +2.75 | +2.17 | +2.46 | -3.79 | -3.53 |

The table shows one consistent pattern. In control cells, PGE2 and adenosine raise TNF/NF-kB (nuclear factor kappa B), inflammatory and interferon-gamma programmes and lower cell-cycle programmes (E2F targets, G2/M checkpoint). In the Crtc2 line these programmes go the other way: lower inflammatory and interferon signatures and higher cell-cycle signatures than control, already without treatment and under PGE2 and adenosine. Under TNF, the TNF/NF-kB programme itself is lower in the Crtc2 line.

**Top over-represented GO biological-process terms** (three strongest per direction):

| Comparison | Higher in first group | Lower in first group |
| :--- | :---: | :---: |
| Crtc2-un vs control-un | G protein-coupled receptor signaling pathway; leukocyte migration; mononuclear cell migration | G protein-coupled receptor signaling pathway; regulation of immune effector process; positive regulation of cell adhesion |
| Crtc2-TNF vs control-TNF | monocarboxylic acid metabolic process; fatty acid metabolic process; lipid biosynthetic process | cytokine-mediated signaling pathway; response to bacterium; regulation of immune effector process |
| Crtc2-PGE2 vs control-PGE2 | chromosome segregation; nuclear division; nuclear chromosome segregation | immune effector process; regulation of immune effector process; response to bacterium |
| Crtc2-adenosine vs control-adenosine | chromosome segregation; nuclear division; nuclear chromosome segregation | adaptive immune response; regulation of immune effector process; immune effector process |
| control-PGE2 vs control-un | immune effector process; negative regulation of immune system process; cytokine-mediated signaling pathway | chromosome segregation; nuclear division; nuclear chromosome segregation |
| control-adenosine vs control-un | T cell activation; cytokine-mediated signaling pathway; defense response to protozoan | chromosome segregation; nuclear chromosome segregation; sister chromatid segregation |

**Number of significant terms (FDR < 0.05)**

| Comparison | GO BP up | GO BP down | KEGG up | KEGG down | GSEA GO BP | GSEA KEGG | GSEA Hallmark |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| Crtc2-un vs control-un | 113 | 533 | 6 | 17 | 266 | 10 | 15 |
| Crtc2-TNF vs control-TNF | 137 | 516 | 14 | 18 | 75 | 24 | 14 |
| Crtc2-PGE2 vs control-PGE2 | 177 | 1,067 | 12 | 38 | 479 | 32 | 18 |
| Crtc2-adenosine vs control-adenosine | 186 | 604 | 26 | 28 | 465 | 4 | 18 |
| control-PGE2 vs control-un | 1,024 | 253 | 43 | 12 | 1,190 | 102 | 28 |
| control-adenosine vs control-un | 273 | 203 | 20 | 4 | 1,008 | 100 | 25 |

Complete tables and dot plots are in `Enrichment/<comparison>/{GO,KEGG,GSEA}/`.

### 7.4 How the Crtc2 line responds to each treatment

The client's comparisons test the Crtc2 line against control within each treatment. To show directly how much of each treatment response remains in the Crtc2 line, we took the genes that respond in control cells and compared the size of their change (treated vs untreated, mean of three replicates) in the Crtc2 line with the size in control. The fraction below is the slope of that relationship: 100% would mean the Crtc2 line responds exactly like control, 0% that it does not respond at all. These are descriptive summaries, not a statistical test (see Section 8.3).

| Treatment | Genes used | Crtc2-line response as a fraction of control | Correlation | Genes responding less than half as much |
| :--- | :---: | :---: | :---: | :---: |
| PGE2 | 2,174 | 58% | 0.86 | 33% |
| adenosine | 254 | 32% | 0.46 | 55% |
| TNF | 3,950 | 74% | 0.92 | 33% |

Genes used: for PGE2 and adenosine, the DEGs of comparisons 5 and 6. TNF vs untreated was not among the requested comparisons and was not tested, so for TNF the genes are those whose mean expression differs at least 2-fold between control-TNF and control-un. Values: `DE_PCA_Results/Crtc2_line_treatment_response_retention.csv`.

The Crtc2-line differences also point in a consistent direction relative to the treatments. Without treatment, of the 740 genes that differ both between the cell lines and with PGE2 in control cells, 91% change in opposite directions (genome-wide correlation of the two sets of fold changes -0.45); for adenosine the corresponding figures are 90% of 146 genes (correlation -0.09). Under PGE2, 97% of the 1,384 genes significant in both comparison 3 and comparison 5 move in opposite directions (correlation -0.77): the difference between the lines under PGE2 is, for the most part, the PGE2 response that the Crtc2 line does not make.

**Selected genes** (log2FC, first group relative to the second; * = statistically significant and at least 1.2-fold; n.t. = not tested, expression too low). Nr4a2, Crem, Pde4b, Pde4d and Sik1 are classical CREB/CRTC target genes; Il6, Ptgs2, Cxcl2, Cxcl10 and Ccl5 are inflammatory genes.

| Gene | Crtc2-un vs control-un | Crtc2-TNF vs control-TNF | Crtc2-PGE2 vs control-PGE2 | Crtc2-adenosine vs control-adenosine | control-PGE2 vs control-un | control-adenosine vs control-un |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: |
| Crtc2 | -0.58* | -0.93* | -0.60* | -0.59* | +0.38* | +0.01 |
| Nr4a2 | -0.40* | -1.04* | -0.85* | -0.80* | +0.85* | -0.03 |
| Crem | -0.67* | -1.42* | -0.57* | -0.73* | +0.19 | -0.10 |
| Pde4b | -0.59* | -1.88* | -1.29* | -0.68* | +1.00* | -0.00 |
| Pde4d | -1.19* | -2.55* | -2.22* | -1.87* | +2.17* | +0.08 |
| Sik1 | -0.66* | -1.08* | -0.91* | -0.60* | +0.37* | -0.01 |
| Il6 | -0.24 | -3.37* | -2.37* | -0.54* | +2.69* | +0.01 |
| Ptgs2 | n.t. | -2.48* | -2.06* | n.t. | +2.33* | n.t. |
| Cxcl2 | n.t. | -0.63* | -0.98* | n.t. | +2.14* | n.t. |
| Cxcl10 | -0.70* | -1.54* | -1.46* | -0.19 | +0.58* | -0.01 |
| Ccl5 | -0.87* | -1.52* | -1.83* | -1.43* | +2.24* | +0.95* |
| Nfkbia | +0.10 | -1.17* | -0.39* | -0.05 | +0.67* | +0.27* |
| Tnfaip3 | +0.02 | -0.76* | -0.41* | -0.10 | +0.48* | +0.04 |
| Fos | +1.59* | +2.34* | +0.16 | +1.56* | +0.19 | +0.06 |
| Fosb | +1.35* | +3.00* | +0.39* | +1.75* | +0.05 | +0.01 |

The CREB targets rise with PGE2 in control cells and are lower in the Crtc2 line in every condition, most strongly under PGE2 and TNF. Il6, Ptgs2 and Cxcl2 are induced by PGE2 in control cells and markedly lower in the Crtc2 line under PGE2 and TNF. Fos and Fosb go the other way: they are higher in the Crtc2 line without treatment and under TNF and adenosine. Values: `DE_PCA_Results/CREB_CRTC_target_and_inflammatory_marker_panel.csv`.

## 8. Conclusions, Limitations and Next Steps

### 8.1 Per-comparison summary

| Comparison | DEGs | Read-out |
| :--- | :---: | :---: |
| Crtc2-un vs control-un | 1,151 | Crtc2 line differs from control already without treatment; CREB targets and inflammatory/interferon genes lower, cell-cycle genes higher |
| Crtc2-TNF vs control-TNF | 5,905 | Largest difference; the Crtc2 line mounts a smaller TNF response, with lower TNF-induced inflammatory genes |
| Crtc2-PGE2 vs control-PGE2 | 2,105 | Mostly the PGE2 response that the Crtc2 line does not make: 66% of these DEGs are PGE2-responsive in control cells, 97% of them moved in the opposite direction |
| Crtc2-adenosine vs control-adenosine | 1,304 | Mostly the baseline Crtc2-line difference: 68% of these DEGs are also DEGs without treatment, all in the same direction |
| control-PGE2 vs control-un | 2,174 | Strong PGE2 response: CREB targets and inflammatory genes up, cell-cycle genes down |
| control-adenosine vs control-un | 254 | Mild adenosine response, same direction as PGE2 but much smaller |

### 8.2 Practical reading

The seven points of Section 2, each with the section and delivered file that supports it.

- **The data are of good quality and the sample labels are consistent.** All 24 libraries passed quality control (87.4–88.8% of reads uniquely mapped), all were sequenced together on one flowcell lane, and replicates of the same condition are 33.2 times closer to each other than to other conditions. No sample needs to be excluded. *(Evidence: Sections 3 and 7.1; `QC/multiqc/`, `DE_PCA_Results/PCA_AllSamples_Overview.pdf`)*

- **The Crtc2 line has 33–47% less Crtc2 mRNA than control, in every condition** (Crtc2 line / control = 0.53–0.67). The transcript is reduced, not absent. RNA-seq measures mRNA only, so whether CRTC2 protein is lost should be confirmed by Western blot (Section 8.4). *(Evidence: Section 7.2 Crtc2 table; `DE_PCA_Results/<comparison>/*_ALL_genes.csv`, row Crtc2)*

- **Without treatment, the Crtc2 line already differs from control in 1,151 genes, and the difference runs opposite to the PGE2 response.** Classical CREB target genes (Nr4a2, Crem, Pde4b, Pde4d, Sik1) are lower in the Crtc2 line; of the 740 genes changed both by the Crtc2 line and by PGE2, 91% change in opposite directions. *(Evidence: Sections 7.3 and 7.4; `DE_PCA_Results/Crtc2-un_vs_control-un/`, `CREB_CRTC_target_and_inflammatory_marker_panel.csv`)*

- **PGE2 is a strong stimulus in control cells (2,174 genes change), and the Crtc2 line reproduces only about 58% of that response.** In control cells PGE2 raises CREB targets and inflammatory genes (Pde4d +2.2, Il6 +2.7 log2 fold change) and lowers cell-cycle genes. Under PGE2 the Crtc2 line differs from control in 2,105 genes; 66% of them are genes that PGE2 changes in control cells, and 97% of those move in the opposite direction. *(Evidence: Sections 7.1 and 7.4; `DE_PCA_Results/control-PGE2_vs_control-un/`, `Crtc2_line_treatment_response_retention.csv`)*

- **Adenosine has a mild effect on control cells (254 genes, only 13 of them changing 2-fold or more) and almost none on the Crtc2 line.** In the Crtc2 line, adenosine-treated and untreated samples are practically indistinguishable (Section 7.1); on the control-cell adenosine genes the Crtc2 line reproduces about 32% of the response. *(Evidence: Section 7.1 (centroid distances) and 7.4; `DE_PCA_Results/control-adenosine_vs_control-un/`)*

- **TNF is the strongest stimulus in the experiment, and the largest Crtc2-line difference is seen under TNF (5,905 genes).** The Crtc2 line still responds to TNF, but less (about 74% of the control response), and TNF-induced inflammatory genes such as Il6 (-3.4 log2 fold change) and Cxcl10 are clearly lower than in control. *(Evidence: Sections 7.1 and 7.4; `DE_PCA_Results/Crtc2-TNF_vs_control-TNF/`, `Enrichment/Crtc2-TNF_vs_control-TNF/GSEA/`)*

- **What we suggest next.** Confirm CRTC2 protein loss by Western blot, and confirm a small panel of the key genes (Pde4d, Nr4a2, Il6) by qPCR (quantitative PCR) in independent cultures. *(Evidence: Section 8.4)*

### 8.3 Limitations

- **The nature of the Crtc2 modification is not stated.** The data show reduced Crtc2 mRNA; they cannot show whether CRTC2 protein is absent or reduced.
- **The response fractions in Section 7.4 are descriptive.** A formal genotype x treatment interaction test was not part of the requested comparisons (the factorial-design box on the form was not ticked). The figures summarise group means and are not p-values.
- **TNF vs untreated was not a requested comparison**, so no TNF-response DEG list is delivered; the TNF statements rest on the PCA and on comparison 2.
- **Three replicates per condition.** Consistent changes are detected well; small or variable effects can be missed, so a gene absent from a DEG list is not shown to be unchanged.
- **Some chemokines are barely expressed in these cells** (for example Cxcl1 and Ccl2 have fewer than 20 reads in every sample) and could not be tested.
- **RNA level only.** Transcript changes do not establish protein or functional changes; pathway terms name a theme, not a demonstrated mechanism.

### 8.4 Suggested confirmation

| Open question | Experiment that would settle it |
| :--- | :---: |
| Is CRTC2 protein lost in the Crtc2 line? | Western blot for CRTC2 in control and Crtc2 cells |
| Are the key gene changes reproducible? | qPCR of Pde4d, Nr4a2, Il6 (and Crtc2) in independent cultures, untreated and with PGE2 and TNF |
| Are the pathway differences functional? | A functional read-out for the main themes, for example a proliferation assay and cytokine secretion (for example IL-6 by enzyme-linked immunosorbent assay, ELISA) |

## 9. Deliverable Files

| Path | Content |
| :--- | :---: |
| `Reads/All_sample_gene_counts.tsv` | Gene-level estimated counts (Salmon), 24 samples — the matrix for re-analysis |
| `Reads/All_sample_gene_tpm.tsv` | Gene-level TPM, 24 samples |
| `Reads/library_metrics_per_sample.csv` | Read pairs, reads assigned to genes, genes detected |
| `QC/multiqc/` | MultiQC quality-control report |
| `DE_PCA_Results/PCA_AllSamples_Overview.pdf` | PCA of all 24 samples |
| `DE_PCA_Results/PCA_coordinates_overview.csv` | PCA coordinates |
| `DE_PCA_Results/SampleDistance_Heatmap.pdf` | Sample-to-sample distances |
| `DE_PCA_Results/normalized_counts_all_samples.csv` | DESeq2-normalised counts |
| `DE_PCA_Results/DEG_summary_all_comparisons.csv` | DEG counts, 1.2-fold and 2-fold cut-offs |
| `DE_PCA_Results/<comparison>/` | Full and significant result tables, volcano plot, top-50 heatmap (6 comparisons) |
| `DE_PCA_Results/Crtc2_line_treatment_response_retention.csv` | Section 7.4 response fractions |
| `DE_PCA_Results/CREB_CRTC_target_and_inflammatory_marker_panel.csv` | Section 7.4 gene panel |
| `Enrichment/<comparison>/GO/` and `Enrichment/<comparison>/KEGG/` | Over-representation, split by direction |
| `Enrichment/<comparison>/GSEA/` | GSEA: GO BP, KEGG, MSigDB Hallmark |
| `Enrichment/Enrichment_summary_all_comparisons.csv` | Number of significant terms per analysis |
| `mouse_Gene_annotation_20260628.xlsx` | Gene annotation reference |

---

*Zhen Gao, PhD, Principal Bioinformatics Scientist, Athenomics*

