# Dual Cut Project — Whole-Genome Sequencing of Ten RPE-1 Clones

**Report Date:** 2026-11-26
**Prepared by:** Zhen Gao, PhD, Principal Bioinformatics Scientist, Athenomics
**Analysis Platform:** Linux HPC server
**Species:** *Homo sapiens* (GRCh38)
**Tissue / Cell:** hTERT RPE-1 cell line, single-cell–derived clones (as stated by the client)

---

## 1. Objectives

This round answers the three questions set out in your questionnaire:

1. **Genome-wide off-target analysis** — did the dual guide RNA (gRNA) pair cut anywhere other than the intended site?
2. **Confirmation of the intended ~200 bp deletion** at the PPP1R2 target, clone by clone.
3. **Any further changes attributable to the CRISPR system** — whether the edited clones carry more new variation than clones that went through the same procedure without acquiring the deletion.

Ten clones were sequenced: three controls (Ctr_1 to Ctr_3) and seven edited clones (Sample_1 to Sample_7).
Abbreviations used throughout: **WGS** = whole-genome sequencing; **gRNA** = guide RNA; **PAM** = protospacer-adjacent motif; **RNP** = ribonucleoprotein (Cas9 protein pre-complexed with gRNA); **SNV** = single-nucleotide variant; **indel** = insertion or deletion; **SV** = structural variant; **PCR** = polymerase chain reaction; **BAF** = B-allele fraction (the share of sequencing reads carrying the non-reference base at a heterozygous site).

---

## 2. Core Conclusions

This section gives the practical answer for each clone in plain language, without the reasoning.
Section 3 lists what we measured, grouped by topic rather than by clone.
Section 8.1 repeats each conclusion with the basis it rests on, Section 8.2 names the delivered file behind
each finding, Section 8.3 states what these data cannot settle, and Section 8.4 pairs each open question
with the experiment that would close it.

| Clone | What you have | What we suggest |
| :--- | :--- | :---: |
| **Ctr_1, Ctr_2, Ctr_3** | No trace of the 202 bp deletion. Both copies of PPP1R2 are intact. | Usable as your comparison clones, with the caveat in Section 8.3 that they were CRISPR-treated as well. |
| **Sample_1** | The deletion is present, but only in about **29%** of the cells in the vial — this is a mixed population, not a pure clone. | Re-clone from a single cell before using it for phenotype work. |
| **Sample_2** | Same situation: the deletion is present in about **32%** of cells. | Re-clone before use. |
| **Sample_3** | One of the two copies carries the deletion; the clone behaves as a clean heterozygote. | Ready to use as a single-allele clone. |
| **Sample_4** | One of the two copies carries the deletion; clean heterozygote. | Ready to use as a single-allele clone. |
| **Sample_5** | **Not two deleted copies.** One copy carries the 202 bp deletion; the other carries a small insertion/deletion at the first cut site. Both copies are disrupted, but by two different events. | Usable if "both copies disrupted" is what you need; not usable if the experiment requires two identical deletions. |
| **Sample_6** | The target region is contaminated by what looks like a PCR product (amplicon) of the same locus: ~16× more reads there than in any other clone. **The genotype cannot be read at face value**, and the contamination also inflates the vector signal for this clone. | Re-extract DNA from a fresh culture and re-check this clone before using it. |
| **Sample_7** | Both copies carry the 202 bp deletion, cleanly. | Ready to use as a two-allele clone. |

**Across all ten clones:** we found **no evidence of off-target cutting** at any of the 1,584 candidate
sites the two gRNAs could plausibly reach, and **no increase in new variation** in the edited clones
relative to the controls.

**One finding that affects how every group comparison in this report should be read:** the seven edited
clones and the three control clones are **not matched in karyotype**. The short arm of chromosome 12
(chr12p) is present in **four copies in all seven edited clones and in two copies in all three controls**.
The most likely explanation is that the two sets were drawn from different subpopulations of the parental
culture. Consequently, any difference between the two groups travels with lineage as well as with
editing, and is reported throughout as **"segregating with the edit"** rather than as "caused by the edit."

---

## 3. Key Findings

1. **The intended deletion is exactly 202 bp and identical in every clone that carries it.** It removes
   chr3:195,542,948–195,543,149, with a single shared base (C) at the join — the signature of direct
   end-joining repair.
2. **No off-target cutting was detected.** Of 1,584 candidate sites, 1,556 showed no evidence in any
   clone; 20 showed the same low-level signal in both edited and control clones (i.e. background); 7 were
   seen only in edited clones and 1 only in a control clone. None reached statistical significance
   (smallest p = 0.18). Read-level inspection of the 7 edited-only sites found **no instance of the same
   indel recurring in more than one clone**, which is what a genuine cut site would produce.
3. **The two gRNAs can also pair up at five PPP1R2 pseudogenes elsewhere in the genome. None of them was cut.**
   All ten clones showed zero junction-spanning reads at those five sites, with 8–26 intact reads per cut
   site confirming the region was well covered.
4. **No new vector integration was detected.** The two integration sites your parental line already
   carries were found in all ten clones, as expected. The plasmid backbone (the AmpR ampicillin-resistance gene) was absent from all
   ten clones, and the Cas9 vector was likewise absent — consistent with your statement that Cas9 was
   delivered as RNP this round.
5. **Edited clones do not carry more new variation than controls** — not for SNVs (p = 0.15), indels
   (p = 0.64), the two combined (p = 0.45), or structural variants (p = 0.38).

---

## 4. Sample Information

| Clone | Role | Expected genotype (your gel/Sanger) | Read pairs sequenced | Raw depth | Usable depth |
| :--- | :---: | :---: | :---: | :---: | :---: |
| Ctr_1 | Control | no deletion | 576,728,998 | 53.8× | 27.6× |
| Ctr_2 | Control | no deletion | 544,862,727 | 50.8× | 23.4× |
| Ctr_3 | Control | no deletion | 497,440,989 | 46.4× | 25.6× |
| Sample_1 | Edited | single-allele deletion | 567,126,162 | 52.9× | 25.6× |
| Sample_2 | Edited | single-allele deletion | 560,871,311 | 52.3× | 22.7× |
| Sample_3 | Edited | single-allele deletion | 596,138,972 | 55.6× | 32.8× |
| Sample_4 | Edited | single-allele deletion | 535,487,760 | 50.0× | 24.6× |
| Sample_5 | Edited | both alleles deleted | 499,253,985 | 46.6× | 24.7× |
| Sample_6 | Edited | both alleles deleted | 675,460,326 | 63.0× | 20.7× |
| Sample_7 | Edited | both alleles deleted | 569,098,815 | 53.1× | 28.2× |

**Two depth figures are given because they answer different questions.** *Raw depth* is total sequence
produced divided by genome size — this is the figure the 45× order refers to, and every clone exceeded it.
*Usable depth* is what remains after removing duplicate reads and after counting each DNA fragment once
rather than twice where the two reads of a pair overlap. The gap between the two is fully explained by the
library insert size (mean 145–187 bp, against a 150 bp read length), which makes the two reads of most
pairs overlap almost completely. All detection statements in this report are based on usable depth
(Figure 1).

Identity checking across 26,674 common variant sites confirmed that all ten samples come from the same
individual (pairwise agreement 99.73–99.81%) and that none carries detectable foreign human DNA
(0.12–0.26% mismatching reads, the level expected from sequencing error alone) — Figure 3.

---

## 5. Analysis Rationale and Decision Criteria

| Question | How it was decided | Threshold used |
| :--- | :--- | :---: |
| Is the deletion present? | Reads that span the join between the two cut sites, counted separately at each cut site. Combining the two cut sites would hide a clone carrying a deletion on one copy and a small indel on the other — which is exactly what Sample_5 turned out to be. | ≥20 bp of the read matching the far side of the deletion, ≥90% sequence agreement |
| How many copies carry it? | Two independent estimates that fail in different ways: sequencing depth across the deleted interval, and the proportion of junction reads among all reads crossing the cut site. | Agreement between the two estimates |
| Was anything cut off-target? | Every site in the genome matching either gRNA with up to 4 mismatches, then read-level inspection at each one. | ≥3 distinct reads showing an indel or a cut-site break, and ≥10% of reads at that position |
| Is a vector integration new? | The genomic breakpoints found in a clone are compared with those shared by all ten clones (inherited from the parental line). | ≥2 independent reads; read support compared against the two known parental sites |
| Do edited clones carry more new variation? | Variants private to one clone — present in that clone and confirmed absent, with adequate coverage, in all nine others. | Quality ≥20, depth ≥10, allele balance ≥0.25; the other nine clones must each have ≥10 reads supporting the reference |

A note on the last row: because each clone was expanded separately from a single cell, every clone
accumulates its own private variants during culture. Attributing any individual variant to CRISPR is not
possible. What can be tested is whether the edited clones, as a group, carry **more** such variants than
the controls. With seven versus three clones the smallest achievable p-value is 0.008, so the test is able
to detect only a complete separation between the groups.

---

## 6. Methods

Reads were aligned to the GATK GRCh38 reference with bwa-mem2 2.2.1 and duplicate-marked with
GATK/Picard 3.4.0 MarkDuplicates. Coverage was computed with mosdepth 0.3.14 in its default mode, which
corrects for the overlap between the two reads of a pair. Adapter trimming was deliberately not performed:
the aligner soft-clips adapter sequence rather than aligning it, which can only cause a true event to be
missed, never create a false one.

Germline SNVs and indels were called jointly across all ten clones with GATK 4.6.2.0 HaplotypeCaller in
GVCF mode followed by GenotypeGVCFs, and filtered with the GATK hard-filter thresholds. Base quality score
recalibration and variant quality score recalibration were not applied, as the required training resources
are not available on this platform; this is noted in Section 8.3.

Structural variants were called with DELLY 2.3.0, first per clone and then re-genotyped across all ten
clones at a unified set of 19,174 sites.

Off-target candidates were enumerated for three spacer versions (gRNA1 at 20 nt, gRNA2 as written at 19 nt,
and gRNA2 extended to 20 nt from the genome) allowing up to 4 mismatches, with the PAM read from the
reference afterwards. The intended target site was required to be recovered by this procedure with zero
mismatches as a positive control before any site was interpreted.

Vector screening used the two plasmid maps from the earlier round. Regions of the reporter vector that
share sequence with the human genome were masked before computing coverage; without this mask, the 53BP1
sequence carried by the vector produces an apparent integration signal at the endogenous 53BP1 locus.

**Positive and negative controls built into the analysis.** The intended target site had to be detected by
the off-target pipeline; the two vector integration sites inherited from the parental line had to be
detected in all ten clones; and the plasmid backbone had to be absent. All three held.

---

## 7. Results

### 7.1 The intended deletion

The deletion removes **chr3:195,542,948–195,543,149 (202 bp)**, joining at chr3:195,542,947 and
chr3:195,543,150 with one shared base (C) at the junction. The same breakpoints were found in every clone
carrying the deletion.

| Clone | Junction reads | Deleted-allele fraction (cut 1 / cut 2) | Reading |
| :--- | :---: | :---: | :---: |
| Ctr_1 / Ctr_2 / Ctr_3 | 0 / 0 / 0 | 0.00 / 0.00 | No deletion; 12–25 intact reads per cut site confirm coverage |
| Sample_1 | 6 | 0.29 / 0.32 | Deletion present in ~29% of cells |
| Sample_2 | 4 | 0.33 / 0.24 | Deletion present in ~32% of cells |
| Sample_3 | 20 | 0.54 / 0.61 | One of two copies |
| Sample_4 | 8 | 0.67 / 0.50 | One of two copies |
| Sample_5 | 13 | 1.00 / 0.50 | One copy deleted; the other carries a small indel at cut 1 (12 reads) |
| Sample_6 | 8,379 | 1.00 / 1.00 | Not interpretable — see below |
| Sample_7 | 22 | 1.00 / 1.00 | Both copies deleted |

**Sample_6.** The target window contains 19,620 reads against 982–1,290 in every other clone, and
sequencing depth *rises* inside the supposedly deleted interval instead of falling. This is the signature
of a PCR product of this locus contaminating the DNA preparation, not of a genomic genotype. Notably, the
same contamination raises this clone's apparent vector copy number (Section 7.4). Its deletion status
should be re-assessed on freshly extracted DNA. All other clones showed 761–997 reads in the same window,
i.e. no detectable cross-contamination.

**Sample_1 and Sample_2** were each judged sub-clonal by two estimates that detect the deletion with very
different efficiency — depth across the interval, and the share of junction-spanning reads. Both gave the
same answer (0.29 and 0.32), which is why the result is not attributed to a quirk of either method.

### 7.2 Off-target screen

1,587 candidate sites were enumerated across the genome; 1,584 are off-target (the remaining 3 are the
intended site, recovered by all three spacer versions as required). Median coverage at candidate sites was
28 reads.

| Evidence pattern | Sites |
| :--- | :---: |
| No evidence in any clone | 1,556 |
| Present in both edited and control clones (background) | 20 |
| Edited clones only | 7 |
| Control clones only | 1 |

The smallest p-value across all 1,584 sites was 0.18; with seven edited and three control clones the
smallest value obtainable is 0.008, so nothing approached significance. The seven edited-only sites were
examined read by read (395 supporting reads in total): **in no case did the same indel appear in more than
one clone**, which is the pattern a real cut site produces. All seven carry four mismatches to the guide —
the weakest class in the candidate list. We therefore read them as clone-specific background rather than
as off-target activity (Figure 5).

Because your experiment uses two guides, we also asked whether the pair could act together anywhere else.
Across the genome there are only **six** places where a gRNA1 and a gRNA2 match site lie 50–500 bp apart:
the intended target and five PPP1R2 pseudogenes (on chromosomes 5, 6, 13, 13 and 21). **All five
pseudogenes showed zero junction reads in all ten clones**, with 8–26 intact reads spanning each cut site
— a negative result with coverage behind it, not an absence of data.

### 7.3 Structural variants

Joint structural-variant calling produced 19,174 sites across the ten clones. Of these, **118 were private
to a single clone**: 6, 28, 19, 14, 4, 17 and 6 in the edited clones against 9, 13 and 2 in the controls
(p = 0.38). Two of the 118 — both breakpoints of one inversion in Sample_4 on chromosome 19 — fall within
600 bp of an off-target candidate. We do not read this as off-target activity: a Cas9 cut produces
breakpoints within about 10 bp of the cut site, not 550–590 bp away, and the surrounding region showed the
same read behaviour in all ten clones (237–296 reads, 0–14 split reads, with no separation between edited
and control clones).

### 7.4 Vector sequence

| Measurement | Result | Reading |
| :--- | :--- | :---: |
| Two integration sites inherited from the parental line | Found in all 10 clones | Positive control passed; these are pre-existing, not new |
| Plasmid backbone (AmpR) in the reporter vector | 0.00 copies in all 10 clones | Negative control passed; no backbone carried in |
| Cas9-EGFP vector | 0.02–0.43 copies (Sample_6: 2.71) | Consistent with your statement that Cas9 was supplied as RNP |
| Reporter vector, whole-vector average | 0.28–0.80 copies (Sample_6: 7.09) | Background level; Sample_6 inflated by the contamination in Section 7.1 |

143 vector–genome junctions were detected in total, 13 of them shared by all ten clones. Junctions seen in
only one clone are supported by a median of 2 reads, against a median of 16 reads for the two known
parental sites, so they sit at the level where sequencing artefacts occur. Seven single-clone junctions
have 5 or more supporting reads and are listed in the delivered table for completeness; four of those seven
are in Sample_6 (Figure 6).

### 7.5 Copy number and new variants

Across the genome, one region differs systematically between the two groups: **chr12p is present at four
copies in all seven edited clones and two copies in all three controls.** This was established twice, by
methods that fail in different ways: sequencing depth relative to the ten-clone median, and allele ratios
at heterozygous sites. The allele-ratio evidence is shown in Figure 4 — median BAF 0.232–0.260 in edited
clones (the 3:1 ratio expected from four copies) against 0.417–0.423 in controls, while chr12q and
chromosome 7 show no difference between the groups.

Joint variant calling produced 7,308,235 variant sites, of which 5,679,637 passed filtering and 45,247
carry a ClinVar annotation. After requiring a variant to be present in exactly one clone and confirmed
absent with adequate coverage in the other nine, **27,644 private variants** remained.

| Group | Private SNVs | Private indels | Total |
| :--- | :---: | :---: | :---: |
| Controls (n = 3) | 622–1,463 | 644–2,353 | 1,266–3,816 |
| Edited (n = 7) | 836–2,129 | 707–3,403 | 1,795–5,532 |

Neither comparison is significant (SNVs p = 0.15, indels p = 0.64, total p = 0.45). The highest counts in
both groups belong to the clones with the highest duplicate rates (Ctr_2 at 12.1%, Sample_2 at 15.2%,
Sample_6 at 34.7%), so part of the spread between clones reflects library preparation rather than biology
(Figure 7).

---

## 8. Conclusions, Limitations and Next Steps

### 8.1 Conclusions by clone, with the basis for each

| Clone | Conclusion | Basis |
| :--- | :--- | :---: |
| Ctr_1–Ctr_3 | No deletion on either copy | 0 junction reads with 12–25 intact reads per cut site |
| Sample_1 | Single-allele deletion, ~29% of cells | Two independent estimates agree (0.29 and 0.29–0.32) |
| Sample_2 | Single-allele deletion, ~32% of cells | Two independent estimates agree (0.32 and 0.24–0.33) |
| Sample_3 | Single-allele deletion, clean | Fractions 0.54/0.61, 20 junction reads |
| Sample_4 | Single-allele deletion, clean | Fractions 0.67/0.50, 8 junction reads |
| Sample_5 | One copy deleted + one copy with a small indel at cut 1 | Cut 1: 0 intact reads, 12 reads carrying an indel; cut 2: 13 intact reads |
| Sample_6 | Genotype not determinable from this DNA preparation | 19,620 reads in the window; depth rises inside the deleted interval |
| Sample_7 | Both copies deleted | Fractions 1.00/1.00, 22 junction reads, 0 intact reads at either cut site |
| All ten | No off-target cutting detected | 1,584 candidates, smallest p = 0.18, no recurrent indel among the 7 edited-only sites |
| All ten | No new vector integration detected | Backbone 0.00 copies; the 2 inherited sites present in all 10 |
| Edited vs control | No increase in new variation | SNVs p = 0.15, indels p = 0.64, SVs p = 0.38 |

### 8.2 Practical reading, with the delivered file behind each statement

| Statement | Where the evidence is |
| :--- | :--- |
| The deletion is 202 bp and identical across clones | `on_target/junction.tsv`, `on_target/ontarget_summary.tsv` |
| Sample_1 and Sample_2 are mixed populations | `on_target/ontarget_summary.tsv`, columns `deleted_allele_frac_CUT1/2` |
| Sample_5 carries two different events, not two deletions | `on_target/indels_near_cutsites.tsv` |
| Sample_6 is contaminated at the target locus | `on_target/ontarget_summary.tsv`, column `reads_in_window` |
| No off-target cutting | `off_target/offtarget_evidence_n10.tsv`, `off_target/treated_only_sites_reads.tsv` |
| The five pseudogene sites were covered and negative | `off_target/paired_cut_offtargets_n10.tsv` |
| No new vector integration | `vector/stage3_summary.json`, `vector/integration_sharing_n10.tsv` |
| chr12p differs between the groups | `copy_number/chr12p_baf_crosscheck.tsv` |
| No increase in private variants | `variants/private_burden_per_sample.tsv`, `variants/private_sv_per_sample.tsv` |

### 8.3 What these data cannot settle

1. **Usable depth is about half the raw output** (21–33× against 46–63×), because the library inserts are
   close to the read length. Detection statements hold at the usable figure.
2. **Off-target candidates were enumerated by mismatch only.** Sites that the guide could reach by looping
   out a base (a "bulge") are not covered by this screen.
3. **The control clones are not CRISPR-naive.** They received the same treatment and simply show no
   deletion. A difference seen only in edited clones therefore travels with the edit; it is not evidence of
   causation by the edit.
4. **The two groups differ in karyotype** (chr12p, Section 7.5). Group differences carry a lineage
   component that cannot be separated from an editing component in this design.
5. **Structural variants were called by one method.** Breakpoint-level statements are reported as
   candidates rather than as confirmed structures.
6. **Sample_6's target genotype is not readable** from this DNA preparation, and depth-based reasoning
   fails for that clone at that locus.
7. **Sample_1 and Sample_2 are sub-clonal**, so any measurement made on them is an average over a mixed
   population.
8. **Clone mix-ups between clones of the same parent cannot be detected** by the identity check, which
   confirms only that all samples come from the same individual.
9. **Base quality score recalibration was not applied** (no training resources on this platform), so
   variant quality values are as produced by the sequencer.

### 8.4 Suggested next steps

| Open question | Experiment that would answer it |
| :--- | :---: |
| Is Sample_6 edited on both copies? | Re-extract DNA from a fresh culture of Sample_6 and repeat PCR plus Sanger sequencing across the target, keeping the PCR product physically separate from the DNA preparation |
| Are Sample_1 and Sample_2 usable? | Re-clone each from a single cell, then genotype the daughter clones |
| What is the second allele in Sample_5? | Sanger-sequence across cut 1 to read the exact indel |
| Does the chr12p difference matter for your phenotype? | Karyotype or digital PCR on the parental culture and on both sets of clones |
| Would a deeper off-target search change the answer? | A higher sequencing depth on the current clones would raise sensitivity for low-frequency events at the same candidate sites |

---

## 9. Deliverable Files

| Folder | Contents |
| :--- | :--- |
| `Dual_Cut_Project_WGS_1126.md` | This report — the entry point; every other file below is referenced from this section |
| `figures/` | Figures 1–7 (PDF) |
| `qc/` | Per-clone alignment statistics, depth summaries, MultiQC report |
| `identity/` | Pairwise genotype agreement matrix, contamination estimates |
| `on_target/` | Deletion structure, per-clone evidence table, junction reads, indels at cut sites |
| `off_target/` | Candidate site list, per-clone evidence, paired-cut sites, read-level inspection of edited-only sites |
| `vector/` | Vector coverage and copy number, integration junctions and their sharing across clones |
| `copy_number/` | chr12p allele-ratio evidence, copy-number segments |
| `variants/` | Joint variant calls, private-variant tables, private structural variants |

**Figures**

| Figure | Shows |
| :--- | :--- |
| Figure 1 | Sequencing depth under both conventions, insert size, duplicate rate |
| Figure 2 | Deleted-allele fraction and junction reads per clone |
| Figure 3 | Sample identity matrix and contamination estimates |
| Figure 4 | chr12p copy number from allele ratios |
| Figure 5 | Off-target evidence across 1,584 candidate sites |
| Figure 6 | Vector copy number for both plasmids |
| Figure 7 | Private SNVs, indels and structural variants per clone |

---

*Zhen Gao, PhD*
*Principal Bioinformatics Scientist, Athenomics*
