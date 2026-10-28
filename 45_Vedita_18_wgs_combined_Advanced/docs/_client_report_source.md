# Chr11 Locus Junction-Level Analysis with Read Rescue, and DLD1 Clone Comparison

**Report Date:** @@REPORT_DATE@@\
**Prepared by:** Zhen Gao, PhD, Principal Bioinformatics Scientist, Athenomics\
**Prepared for:** Xin Lin and Vedita Anand Singh, Li lab, Boston Children's Hospital / Harvard Medical School\
**Analysis Platform:** Linux HPC server\
**Species:** Homo sapiens (reference genome hg38 / GRCh38, the same alignments as the 2026-10-09 report)\
**Tissue / Cell type:** DLD1-136 and GM12878 cell lines (populations and clones), as named by the client; no further tissue annotation was supplied\
**Target region:** chr11:114,604,000–114,625,500 (GRCh38), analysed within the window chr11:114,554,000–114,675,500 (target ± 50 kb, as in the earlier report)\
**Samples reported:** S01–S10 (DLD1-136 n = 5, GM12878 n = 5); S11–S18 were used only as a background panel (Section 5)\
**Relation to earlier work:** this report extends Section 7.3 of the report dated 2026-10-09 (`Chr11_Repeat_Region_WGS_SV_Report_1009.md`); that report is unchanged and its findings stand.

---

## 1. Objectives

Your request of 2026-10-13 asked for two analyses on the existing hg38 alignments.

**Analysis 1 — junction-level analysis of the target region with read rescue:** re-extract every read that belongs to the locus, including reads with mapping quality (MAPQ) 0, multi-mapping reads, supplementary and secondary alignments, soft-clipped reads and unmapped reads whose mates map to the region; reassemble locally; and report

1. base-pair breakpoints, including any microhomology or inserted sequence at each junction;
2. junction-spanning read counts per event, before and after read rescue;
3. which junctions are well supported and which are ambiguous because of repeat content;
4. which events are called by both structural-variant (SV) callers, Manta and DELLY, and which by only one.

**Analysis 2 — clone-specific event calling:** events present in the +doxycycline (+Dox) DLD1 clones S04 and S05 and absent from the parental population S01 and the control clones S02 and S03, with the filtering criteria stated.

**Deliverables:** SV VCF (Variant Call Format) files from Manta and DELLY, unfiltered and filtered; single-nucleotide variant (SNV) and small insertion/deletion (indel) VCF files; and an annotated per-event table with a column for the supporting caller(s).

All of these are answered in this report. Analysis 2 is reported twice, at two levels of resolution: at the locus with the read-rescue evidence (Section 7.8), and genome-wide from the caller outputs (Section 7.9). Two further SV callers, GRIDSS and SvABA, were added to the Manta/DELLY comparison at the locus, and junction-spanning primers were designed for confirmation. `tables/client_request_coverage.tsv` lists every question and deliverable in your request with the section and file that answers it.

---

## 2. Core Conclusions

This section gives each conclusion in one sentence, what it means for you, and where the evidence is. Section 3 lists the measurements by topic; Section 8.1 sets out the basis of each conclusion; Section 8.2 explains each result per sample group in more detail; Section 8.3 states what the data cannot show; Section 8.4 pairs each open question with an experiment.

1. **S04 (+Dox clone 3): the broken end is resolved — chr11 at 114,605,216–218 joined to telomeric repeat, and distal of it one chr11q homolog is lost to the end of 11q.** What it means: this clone carries a terminal deletion of 11q whose broken end is capped by telomeric repeat. A junction-specific PCR primer pair (Section 7.11) can confirm it on S04 DNA, with S01–S03 as negative controls. Evidence: Section 7.4, `tables/junction_table_target.tsv` (row J0922), `figures/junction_views/Fig4_J0922.png`, `tables/key_junction_read_evidence.tsv`, `locus_BAM/S04_DLD1_c3_dox.locus_rescue.bam`.
2. **S05 (+Dox clone 4): the same chr11q homolog is lost, but the end is different and more complex — a 376 bp fold-back inversion at the step and a chr21 junction 1.2 kb distal of it.** What it means: at the fold-back the broken end has folded back on itself (it is not an insertion); the chr21 junction carries only about 1 kb of chr11 sequence and does not contradict the distal loss; how the two events relate is open (Section 7.4). The chr21 junction can be tested with the J0786 primer pairs (Section 7.11); the fold-back needs a different design (Section 8.4). Evidence: Sections 7.4 and 7.7, rows J0359 and J0786, `figures/junction_views/Fig4_J0359.png` and `Fig4_J0786.png`, `supporting_analyses/chr11_BAF/retained_homolog_S04_vs_S05.tsv`, `tables/caller_concordance_by_sample.tsv`.
3. **S01–S03 remain a valid baseline: none of them carries any of the three +Dox junctions.** What it means: they can serve as negative controls in junction PCR. S04 belongs to a sub-lineage that none of the control clones represents, so events labelled +Dox-specific in S04 cannot be distinguished from variants already present in that sub-lineage before doxycycline treatment (Section 8.3); S05 is not affected. Evidence: Section 7.4 table (reference fragment counts), `tables/junction_by_sample_long.tsv`.
4. **GM12878 clones: no clone-specific rearrangement at the locus above one fragment.** What it means: no event at this locus in S07–S10 calls for follow-up; the junction pair 13 kb proximal of the target, present in all five GM12878 samples, is a feature of the GM12878 line rather than of any clone. Evidence: Sections 7.5 and 7.8, `tables/junction_table_target.tsv` column `GM12878_clone_specific`, `tables/junction_table_window_context.tsv` rows J1365 and J0533.
5. **Most junction candidates at this locus are single reads inside a 22 kb repeat and should not be followed up individually.** What it means: of 527 candidates, 4 have two or more independent fragments after read rescue, and three of these are the +Dox junctions of points 1–2; follow-up effort is best spent on those. Evidence: Sections 7.1 and 7.5, Figure 5, `tables/junction_table_target.tsv` columns `classification` and `claim_level`.
6. **Genome-wide: 13 candidate structural variants are private to one +Dox clone (8 in S04, 5 in S05); none is shared by the two clones and none lies on chr11.** What it means: no second rearrangement common to both +Dox clones was found. Of the two largest events, the 327 kb chr8 duplication in S04 is confirmed by read depth and allele balance, and the chr4 event in S05 is a real junction but not a deletion; the other 11 are candidates for inspection in the whole-genome BAM files and for PCR (Section 8.4). Evidence: Section 7.9, `genome_wide_DLD1_comparison/dld1_dox_specific_events.tsv`, `supporting_analyses/SV_verification/`.

---

## 3. Key Findings

| # | Finding and measurements | Where |
| :--- | :--- | :---: |
| **1** | **The target region lies inside a 22 kb tandem repeat array.** | |
| 1a | chr11:114,603,768–114,625,660 is an array of a 42 bp unit (`GGTAACCACTGTTACCTGGTGGATAATAAGTATTGCCTCGTG`, about 522 copies), with an internal segmental duplication of 90.2% identity. | Section 7.1; Figure 1 |
| **2** | **Read rescue changed the junction evidence in both directions.** | |
| 2a | Per sample, 18,225–39,040 read pairs were rescued for the window; 9,440–20,554 of them were placed elsewhere in the original alignments and list the window only as an alternative hit. | Section 7.2; Figure 2 |
| 2b | The S04 telomere junction had 0 supporting fragments by the standard count and 4 after rescue. | Section 7.4 |
| 2c | 269 of 527 target-region candidates are better or equally explained by an unrearranged repeat copy once each read competes for its best placement. | Sections 7.3, 7.5 |
| **3** | **Three junctions are specific to the +Dox clones; two of them lie inside the depth step.** | |
| 3a | S04: chr11:114,605,216–218 joined to (TTAGGG)n; 4 fragments (7 read ends); allele fraction 0.29 (95% confidence interval, CI, 0.12–0.55). | Section 7.4; Figure 4 |
| 3b | S05: fold-back inversion chr11:114,605,973 / 114,606,349, 2 fragments; chr11:114,607,192–chr21:20,648,241, 3 fragments. | Section 7.4; Figure 4 |
| 3c | Depth steps: S04 at 114,605,200 (95% CI 114,605,000–114,605,400); S05 at 114,605,800 (95% CI 114,605,600–114,606,000). | Section 7.6; Figure 1 |
| **4** | **Both +Dox clones lost the same one of the two chr11 homologs distal of 114.6 Mb.** | |
| 4a | Distal chr11q (114.7–135.0 Mb): 97.8% (S04) and 98.0% (S05) of sites heterozygous in S01–S03 have become homozygous; at 99.78% of the 8,004 sites homozygous in both clones the same allele is retained. | Section 7.7; Figure 6 |
| 4b | Had the clones lost different homologs, they would retain opposite alleles at nearly every one of these sites; the agreement holds in each of eleven 2 Mb windows (≥ 98.96%). The two losses have different broken ends and occurred in different DLD1 sub-lineages, so they are independent events. Whether this homolog is maternal or paternal cannot be determined without parental genotypes; with two independent clones, the same homolog is expected by chance in half of cases. | Sections 7.7, 8.3 |
| **5** | **Callers rarely agree with the read evidence at this locus.** | |
| 5a | Of 527 target-region junctions, none is called by both Manta and DELLY; the S05 chr21 junction is called by DELLY, GRIDSS and SvABA; 29 are called by GRIDSS alone or together with Manta, mostly single-fragment events. | Section 7.5 |
| **6** | **A widespread sequencing artefact explains a large class of earlier candidates.** | |
| 6a | Reads from short DNA fragments run into the adapter and a poly-G tail; the tail aligns to a poly-G tract at chr2:32,916,3xx. This accounts for 693 of the 1,181 earlier target-region candidates. | Section 7.3 |
| **7** | **Genome-wide DLD1 comparison.** | |
| 7a | Of 12,953 DELLY events re-genotyped in S01–S05, 13 are +Dox-specific (8 in S04, 5 in S05) and none is shared by both clones; the largest copy-number change is a 327 kb single-copy gain on chr8 in S04, confirmed by read depth and allele balance. The chr4 event in S05 called as an 888 kb deletion is a rearrangement junction; the 888 kb interval is not lost. The chr11q loss itself produces no caller event. | Section 7.9 |

---

## 4. Sample Information

Sample names, cell lines, doxycycline status and clone numbers are as supplied by the client. "Read pairs in the rescue set" counts read pairs (identified by their read name); Section 5 defines the sources.

| Sample | Cell line | Treatment | Role in this report | Autosomal depth (×, 2026-10-09 report) | Read pairs in the rescue set |
| :--- | :---: | :---: | :---: | :---: | :---: |
| S01 DLD1_pop | DLD1-136 | no Dox | DLD1 baseline (population) | 16.4 | 25,154 |
| S02 DLD1_c1_ctrl | DLD1-136 | no Dox | DLD1 baseline (control clone 1) | 14.4 | 24,336 |
| S03 DLD1_c2_ctrl | DLD1-136 | no Dox | DLD1 baseline (control clone 2) | 15.9 | 28,550 |
| S04 DLD1_c3_dox | DLD1-136 | +Dox | case (clone 3) | 16.4 | 24,236 |
| S05 DLD1_c4_dox | DLD1-136 | +Dox | case (clone 4) | 18.8 | 33,654 |
| S06 GM12878_pop | GM12878 | no Dox | GM12878 baseline (population) | 12.9 | 18,225 |
| S07 GM12878_c5 | GM12878 | no Dox | clone 5 | 15.9 | 26,505 |
| S08 GM12878_c6 | GM12878 | no Dox | clone 6 | 15.2 | 26,499 |
| S09 GM12878_c7 | GM12878 | no Dox | clone 7 | 18.8 | 34,400 |
| S10 GM12878_c8 | GM12878 | no Dox | clone 8 | 18.5 | 39,040 |

S11–S18 (Mutu I and DG75) were processed identically and used only to flag junctions that also occur in unrelated cell lines; they are not reported sample by sample. The genome-wide comparison (Section 7.9) uses S01–S05; the genome-wide SNV/indel calls (Section 7.10) cover S01–S10.

---

## 5. Analysis Rationale and Decision Criteria

**What is examined.** All alignment records stored within the analysis window chr11:114,554,000–114,675,500, plus every read pair that belongs to the locus but is stored elsewhere. A read pair enters the rescue set when (a) any of its records lies in the window, (b) its mate is placed in the window, (c) a supplementary alignment (the SA tag: another segment of the same read aligned elsewhere) lies in the window, (d) an alternative alignment listed by the aligner (the XA tag: an equally or nearly equally good placement) lies in the window, or (e) both mates are unmapped and one aligns to the window sequence. For every rescued pair, all of its records anywhere in the genome are retrieved. Nothing is selected from the earlier candidate list; earlier events are only compared afterwards.

**Counting units.** A *read end* is one sequenced read; a *fragment* is one read pair (one DNA molecule). Alternative alignments of the same read are kept as positional evidence but never counted as additional support. Reads marked as PCR or optical duplicates are kept in the files and excluded from counts. All thresholds below use fragments.

**Before and after rescue.** *Before* reproduces the earlier standard count: non-duplicate split reads whose two parts both have MAPQ ≥ 20, plus discordant pairs whose mates both have MAPQ ≥ 20. *After* re-tests every rescued read of every sample against a template of the junction sequence and against its best placement anywhere in the genome; a read end supports the junction only if it crosses the junction by at least 20 bases on each side and its alignment score (AS) on the junction template exceeds its best genome alignment by at least 10. A discordant pair supports it when its mates lie on the correct sides within the fragment-size limit (500 bp).

**Evidence class (fixed before the results were examined).**

| Class | Condition | Meaning |
| :--- | :---: | :---: |
| Well supported | ≥ 3 fragments in one sample; ≥ 2 kinds of evidence (split reads, discordant pairs, assembled contig); both 100 bp junction flanks placeable at a single genomic position (no other placement within 10 AS points); breakpoint uncertain by ≤ 10 bp | Ready for PCR confirmation |
| Limited | ≥ 2 fragments; flanks placeable at a single position; uncertainty ≤ 10 bp; but < 3 fragments or one kind of evidence | Real candidate with modest support |
| Ambiguous (repeat) | ≥ 2 fragments, but a flank, or the partner end, can be placed equally well at more than one position, or the breakpoint can shift by > 10 bp within identical sequence | Junction present, location of one end not unique |
| Weak | 1 fragment | Listed, not interpreted |
| Not supported | 0 fragments after rescue | Read evidence is explained at least as well without the junction |

**Claim level (what may be concluded from the evidence).** The evidence class describes the reads; the claim level describes what a reader may conclude. It adds one test that uses information independent of the junction reads themselves, the copy-number step and the B-allele fraction (BAF), but still from the same sequencing data.

**Coherence with the copy-number step.** A junction is called *coherent* in a sample when two tests, both computed in that sample, are met: (1) *position*: one of its chr11 breakpoint intervals overlaps the 95% confidence interval of that sample's depth step; and (2) *loss*: at least 90% of the sites heterozygous in all of S01–S03 are homozygous in that sample distal of the step (114.7–135.0 Mb), i.e. one copy of distal 11q is really lost. Test (2) is needed because a sample without a step has a confidence interval spanning the whole window, so test (1) alone would be met by any junction. Coherence uses only the chr11 end of the junction; it says that the junction lies where one copy is lost, not that it is the only structure there, and it is not independent confirmation, because the junction reads, the depth and the allele fractions all come from the same sequencing run. The tests were applied to every junction with two or more fragments in a sample (`tables/structural_coherence_tests.tsv`).

| Claim level | Condition | Wording used in this report |
| :--- | :---: | :---: |
| Mapping observation | ≤ 1 fragment | not interpreted |
| Repeat-ambiguous junction | Ambiguous (repeat). Exception: if only the partner end is ambiguous and the chr11 end is unique, the junction can still be coherent (definition above) and is then classed in the last-but-one row | junction-like evidence with ambiguous placement |
| Limited junction candidate | Limited, not coherent | limited junction candidate |
| Well-supported junction candidate | Well supported, not coherent | well-supported junction candidate |
| Limited, coherent with the copy-number step | Limited, Well supported, or Ambiguous only at the partner end; and coherent with the copy-number step in a sample with ≥ 2 fragments (both tests of the definition above) | limited short-read evidence coherent with the local copy-number pattern |
| Confirmed | any of the above plus an independent assay (PCR or FISH, fluorescence in situ hybridisation) | confirmed |

No independent assay result was available, so no junction can be "confirmed"; this is a statement about the inputs, not a negative result (`tables/external_input_status.tsv`).

**Clone comparison at the locus.** In a sample, an event is *present* with ≥ 1 fragment (*confident* with ≥ 2); *absent* when no fragment supports it and ≥ 10 fragments span the same position without it (an allele at fraction 0.25 would be missed with probability ≤ 0.75^10 ≈ 6%); otherwise *underpowered*. An event is *+Dox-specific* when present in S04 or S05 and absent in all of S01–S03. No allele-fraction minimum is applied; the allele fraction and its 95% CI are reported per sample.

**Clone comparison genome-wide.** Every DELLY site found in any of S01–S05 was re-genotyped in all five samples. Per sample, *alt* is DELLY's count of variant-supporting read pairs plus junction reads (fields DV + RV), *ref* the corresponding reference count (DR + RR); these are read counts as reported by DELLY, not de-duplicated fragments. *Present*: non-reference genotype and alt ≥ 3. *Low support*: alt 1–2 (reported, never used to call an event present or absent). *Absent*: alt = 0 and ref ≥ 10. *Underpowered*: alt = 0 and ref < 10. *+Dox-specific*: present in S04 or S05, alt = 0 in all of S01–S03 and absent in all three. DELLY's own somatic filter (S04, S05 as cases; S01–S03 as controls) is reported as a separate column; it does not change the labels.

**Artefacts excluded (each counted in `tables/evidence_filter_impact.tsv`).** Soft clips and split alignments that consist of adapter read-through or follow it; segments with very low sequence complexity (base-composition entropy < 1.2 bits, e.g. poly-G tails); alternative alignments scoring more than 10 AS points below the read's best alignment.

---

## 6. Methods

**Input.** The coordinate-sorted, duplicate-marked Binary Alignment Map (BAM) files of the 2026-10-09 report (bwa-mem2 v2.2.1 `mem -Y`, GATK hg38 analysis set `Homo_sapiens_assembly38.fasta`), read only.

**Read rescue.** Each BAM was streamed once with samtools 1.23.1 (`view -e`) to collect records whose mate, SA segment or XA alternative hit falls in the window; window records were taken by region query; both-unmapped pairs were aligned to the window sequence with minimap2 2.30 (`-x sr`; kept when ≥ 40 aligned bases, ≥ 95% identity, entropy ≥ 1.6 bits). All records of every rescued read name were then fetched by position. Every alignment record in the window has full read sequence (no hard clipping; verified for all 18 samples).

**Local re-alignment of rescued reads.** The full sequences of rescued reads were re-aligned to the same reference with bwa-mem2 v2.2.1 `mem -a -Y`, which writes every alternative placement as an explicit secondary alignment. This re-aligns the rescue set only, not the sequencing run.

**Junction discovery.** Split reads, soft clips ≥ 20 bp and discordant pairs (both views: original and re-aligned) with at least one end in the window were clustered across all samples with a 15 bp tolerance. Clusters with ≥ 2 fragments were assembled locally with SPAdes v4.3.0 (`--sc --only-assembler`, k = 21, 33, 55; reads trimmed at the adapter) per sample and per cell line; contigs and soft-clip consensus sequences were aligned with bwa-mem2 `mem -a`.

**Genotyping.** For each junction a 400 + 400 bp junction template and two 801 bp reference templates were built; all rescued reads of all samples were aligned to them (bwa-mem2 `mem -a -Y`) and scored as in Section 5. Junctions whose support comes from the same reads and differ only in the placement of the partner end were merged into one event, listing the alternative placements.

**Breakpoint annotation.** Microhomology and breakpoint uncertainty were computed from the reference by sliding the junction while the junction sequence stays identical; inserted sequence comes from the reads or contigs. For every supporting read end, the overlap (positive, microhomology) or gap (negative, inserted bases) between its consecutive aligned segments was also measured directly from the alignments (`tables/key_junction_read_evidence.tsv`). Each 100 bp flank was aligned to GRCh38 with `mem -a` to count equally good placements. Repeat annotation: UCSC hg38 RepeatMasker, simpleRepeat (tandem repeats detected by Tandem Repeats Finder) and genomicSuperDups (segmental duplications) tables.

**Caller concordance.** Four callers were compared with every junction. Manta 1.6.0 (per sample genome-wide `diploidSV` and `candidateSV`; joint target-window run) and DELLY v2.3.0 (per sample discovery; window re-genotyping) are the outputs of the 2026-10-09 run. GRIDSS (joint run over all 18 samples) and SvABA (joint run over S01–S05) were run on the analysis window with the same BAM files. A caller supports a junction when one of its records has both ends within 100 bp of the junction's two ends (as in the earlier report) with the same orientation; records with only one end nearby, and GRIDSS single-breakend records, are listed as anchors and do not count as support. The smallest end-to-end distance of each supporting record is reported, so a stricter tolerance can be applied from the table.

**Primer design.** For every junction with at least two supporting fragments, junction-spanning primer pairs (product 150–350 bp, melting temperature 57–63 °C) were designed with Primer3 (primer3-py 2.4.0) on the junction template, first with both primers outside the junction and then with one primer required to cross it. Each primer and its 3′-terminal 15 bases were counted as exact matches on both strands of the whole reference, and every full-length match was located (`tables/primer_binding_sites.tsv`). A pair is *strict* when every primer and every 3′ 15-mer occurs at most once, and *junction-specific* when one primer cannot prime the reference at all (no full-length match, at most one 3′ 15-mer match) and the other has exactly one full-length reference match, at the expected flank.

**Depth step and BAF.** Depth per base (samtools `depth`, MAPQ ≥ 20 and ≥ 0, duplicates excluded) in 200 bp bins, normalised to each sample's own median over chr11:114.40–114.58 Mb and divided by the mean of S01–S03; a two-level step model with 1,000 block bootstraps gave the boundary and its 95% CI. BAF = alternative-allele reads / all reads at sites heterozygous in every baseline sample of the cell line (allele fraction 0.25–0.75, depth ≥ 10), from the joint SNV calls of Section 7.10 restricted to chr11.

**Genome-wide SNV/indel calls.** S01–S10 were called jointly with bcftools 1.23.1 (`mpileup -q 20 -Q 20 --max-depth 1000 -a AD,DP`, `call -mv`, `norm -m -any`) on chromosomes 1–22, X, Y and the mitochondrial genome, in 5 Mb blocks. Two files are delivered: all calls, and calls with variant quality (QUAL) ≥ 30. All sites were called as diploid, including chrY and the mitochondrial genome; genotypes there should be read as presence of the alternative allele.

**Genome-wide DLD1 clone comparison.** The per-sample DELLY discovery sites of S01–S05 were merged (`delly merge`, quality and allele-fraction gates relaxed so that the site list is not pre-filtered) and every site was re-genotyped in each of the five BAM files (`delly sr -v`, same exclusion of alternate, decoy and HLA contigs as in the earlier report); DELLY's somatic filter was applied in tagging mode, so no record is removed. Each event was matched to the per-sample Manta calls (both ends within 100 bp, same orientation) and annotated with GENCODE v45 protein-coding and lncRNA genes containing either breakend.

**Delivered file headers.** Local directory paths were removed from the headers of all delivered BAM and VCF files; record counts and a checksum of the records were compared before and after (`methods/header_sanitization_audit.tsv`). Every DELLY VCF was compared with its source BCF in the same way (`methods/delly_conversion_audit.tsv`).

---

## 7. Results

### 7.1 The target region is a tandem repeat array

The simpleRepeat annotation places a 42 bp unit repeated about 522 times across chr11:114,603,768–114,625,660, which covers the whole target region (Figure 1, top). Adjacent copies differ by a few bases, so a 150 bp read usually still has a best position, but often only a few alignment-score points ahead of the next copy. This is why the earlier report flagged many events as repeat-ambiguous, and why junction support has to be judged by competition between placements rather than by mapping quality alone.

![Figure 1 — locus overview](figures/Fig1_locus_overview.png)

*Figure 1. Top: target region and the tandem repeats annotated by simpleRepeat (repeat unit ≥ 20 bp). Orange: the 42 bp repeat family, i.e. the array of about 522 copies (114,603,768–114,625,660), a 47 bp variant immediately proximal of it and a short distal remnant (114,631,682–114,632,022). Grey hatched: an AT-rich low-complexity block (114,601,509–114,603,053; repeat-unit GC ≤ 22%) that is not part of the 42 bp family; depth there is near zero in all 18 samples, a coverage gap common to every cell line. Middle: depth of S04 and S05, and of the no-Dox samples S01 (population, dashed) and S02, S03 (clones), relative to the mean of S01–S03, 200 bp bins, MAPQ ≥ 20. Because S01–S03 form the baseline, their three ratios average 1 in every bin; whether a baseline sample itself has a step is given by its own step fit (`supporting_analyses/depth_boundary/boundary_estimates.tsv`). Shaded bands are the 95% CIs of the S04 and S05 depth steps; the gap at 114.6015–114.603 Mb is the AT-rich block, where the baseline has no MAPQ ≥ 20 reads, and the ratios in the bins at its two edges rest on very few reads. Bottom: every target-region junction with at least one supporting fragment, at the position of its target-side end; colour = class, marker shape = +Dox-specific events. Most points (254 of 258) are Weak junctions with a single fragment, and 247 of these occur in only one sample. They appear at a similar density in every sample, including those without a copy-number step: 21–27 per DLD1 baseline sample (S01–S03), 20–40 per GM12878 sample (S06–S10) and 15–21 in the +Dox clones (S04, S05). This row at the bottom is therefore an internal control: the level of single-fragment, junction-like reads that mis-placement within the repeat produces without any rearrangement. Of the 39 points marked as S04-only or S05-only, 36 are such single-fragment points and are not interpreted; only J0922, J0359 and J0786, with two or more fragments, stand above this background.*

### 7.2 Read rescue

| Sample | Records stored in the window | Read pairs rescued | Of which listed only through an alternative (XA) hit |
| :--- | :---: | :---: | :---: |
| S01 | 26,105 | 25,154 | 12,061 |
| S02 | 24,871 | 24,336 | 11,845 |
| S03 | 29,241 | 28,550 | 13,891 |
| S04 | 22,006 | 24,236 | 13,184 |
| S05 | 30,108 | 33,654 | 18,551 |
| S06 | 17,463 | 18,225 | 9,440 |
| S07 | 26,916 | 26,505 | 12,999 |
| S08 | 26,231 | 26,499 | 13,317 |
| S09 | 34,350 | 34,400 | 17,165 |
| S10 | 36,811 | 39,040 | 20,554 |

About half of each rescue set consists of reads placed elsewhere whose alternative placement lies in the window — the multi-mapping reads you asked us to include. None of the both-unmapped reads (0.25–3.3 million per sample across the 18 samples) aligned to the window sequence. No record in the window of any of the 18 samples is hard-clipped or an explicit secondary alignment, so every record carries its full read sequence (Figure 2).

![Figure 2 — rescue sources](figures/Fig2_rescue_sources.png)

*Figure 2. Read pairs in each sample's rescue set, by the route through which they entered (Section 5): green, pairs with a record stored in the window (almost always with the mate in the window too); purple, pairs placed elsewhere that entered only because an alternative (XA) hit lies in the window. Mates, supplementary segments and unmapped pairs added no further pairs in any sample. Counts are not normalised for sequencing depth, so totals are not comparable between samples.*

### 7.3 An adapter and poly-G artefact

The libraries have a median fragment length of about 160 bp, so many 150 bp reads run through the end of the DNA fragment into the sequencing adapter (`AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC…`) and then into a run of G bases, the signal produced by this chemistry when no template is left. Aligners place that poly-G tail on a 150 bp poly-G tract at chr2:32,916,330–32,916,480 and report the read as split. Across the 18 samples, 16,768 such split alignments, 163,007 adapter soft clips and 668 discordant pairs whose mate is low-complexity sequence were removed from the original alignments (`tables/evidence_filter_impact.tsv`). Of the 1,181 target-region events listed in the earlier report, 693 have a partner end in this tract; the earlier report had already treated it as a mapping artefact, and this analysis identifies the sequence that produces it.

### 7.4 Junctions with two or more fragments in the target region

| Event | Junction (GRCh38) | Microhomology / inserted | Class / claim level | Present in | Fragments before → after rescue | Allele fraction (95% CI) | S01–S03 | Callers |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| J0922 | chr11:114,605,216–218 (+) → (TTAGGG)n | `GT` / none | Ambiguous (partner end) / coherent with step | S04 | 0 → 4 | 0.29 (0.12–0.55) | absent (13, 12, 24 reference fragments) | none |
| J0359 | chr11:114,605,971–973 (+) ↔ chr11:114,606,347–349 (+), 376 bp fold-back inversion | `TA` / none | Limited / coherent with step | S05 | 1 → 2 | 0.17 (0.05–0.45) | absent (15, 11, 28) | none |
| J0786 | chr11:114,607,190–192 (+) → chr21:20,648,241 (−) | `AT` / none | Limited / limited candidate | S05 | 3 → 3 | 0.27 (0.10–0.57) | absent (18, 14, 18) | DELLY, GRIDSS, SvABA |
| J1036 | chr11:114,609,273–279 (−) ↔ 114,609,283 (−), 10 bp inversion | `GTTACC` / none | Limited / limited candidate | S03 | 1 → 2 | 0.09 (0.02–0.27) | S03 only | GRIDSS |

Positions are given as the range within which the junction can slide over identical sequence (the microhomology). "(+)" means the chr11 sequence to the left of the position is retained; "(−)" the sequence to the right. `tables/structural_coherence_tests.tsv` gives the two coherence tests of Section 5 for every sample and event.

**S04 (J0922).** Seven read ends from four DNA fragments read chr11 sequence up to 114,605,216–218 and then only TTAGGG repeats, up to the end of the read or of the DNA fragment (Figure 4, `figures/junction_views/Fig4_J0922.png`). The repeat matches telomeric-repeat tracts at the ends of chrX, chr2 and chr4 and an interstitial tract on chr10 equally well, so the partner end has no unique genomic position; the coordinate printed for it is one representative placement. The chr11 end is unique (one placement for its 100 bp flank). None of the four fragments was counted by the standard method, because the telomeric part of the reads was soft-clipped or aligned elsewhere with MAPQ 0; the junction was found only through the rescued alignments. The junction lies inside the 95% CI of the S04 depth step (Section 7.6). Its allele fraction is compatible with one of two chr11 copies.

**S05 (J0359, J0786).** The fold-back inversion joins the plus strand at 114,605,973 to the minus strand at 114,606,349 (both "+" ends), i.e. the sequence folds back on itself 376 bp further along; it lies inside the 95% CI of the S05 depth step. The chr21 junction (the earlier report's EV0623, refined to the same coordinates and `AT` microhomology; every supporting read end shows the same 2 bp overlap) lies 1.2 kb distal of the step and its chr21 end is uniquely placeable; because it lies outside the step interval it stays at the claim level "limited junction candidate". It is underpowered in S04 (7 reference fragments) because S04 has lost that segment of one copy.

How the chr21 junction relates to the fold-back and to the distal loss is not resolved. Three readings remain open, and the present data do not distinguish them: (i) the segment duplicated by the fold-back broke again and was joined to chr21, i.e. a second round of breakage and rejoining on the same homolog; (ii) the junction is carried by only part of the clone: its allele fraction, 0.27 (3 junction and 8 reference fragments; 95% CI 0.10–0.57), is below the 0.5 expected if every cell carried the junction on one chr11 homolog, but with three junction fragments it cannot be distinguished from 0.5; (iii) it is an event independent of the 11q loss. chr21 copy number was not examined. Junction-spanning PCR with the J0786 primer pairs (Section 7.11) tests whether the junction is real; which homolog carries it and whether every cell carries it would need further assays, for example single-cell PCR.

**S03 (J1036).** A 10 bp inversion inside the repeat array with 2 fragments in the no-Dox control clone S03 only; it is not related to the +Dox phenotype.

![Figure 3 — every junction with two or more fragments](figures/Fig3_key_events_summary.png)

*Figure 3. All 17 junctions in the window (target region ± 50 kb) with two or more fragments in a reported sample (S01–S10), grouped by their pattern across samples: A, +Dox-specific (present in S04 or S05, absent from S01–S03); B, DLD1 no-Dox clone only; C, GM12878 line-wide (present in the GM12878 population S06, absent from DLD1); D, a single GM12878 clone; E, present in the DLD1 baseline and in other cell lines. Within each group, target-region junctions (T) come first, then flank junctions (F), by chr11 position. Cells give independent fragments per sample after rescue; hatched cells have no supporting fragment but fewer than 10 reference fragments, so absence is not established there. The S11–S18 column counts the background samples (Mutu I, DG75) with supporting reads; these samples do not enter the class. Right: claim level and supporting callers.*

![Figure 4 — S04 junction at base resolution](figures/junction_views/Fig4_J0922.png)

*Figure 4 (one of 17; all in `figures/junction_views/`). Junction sequence (wine: chr11 side, dark teal: microhomology or inserted bases, blue: partner side) and every read end that supports it, labelled by sample, cell line, population or clone, fragment and read (R1/R2); "." marks a base identical to the junction sequence. The second header line lists, by cell line, every sample with supporting reads, including the background samples S11–S18 (Mutu I, DG75); the evidence class is computed from S01–S10 only, so a junction that also appears in unrelated cell lines is flagged there.*

### 7.5 All target-region junctions, before and after rescue, and caller concordance

527 distinct junctions have at least one end in the target region (1,847 more lie only in the ± 50 kb flanks; `tables/junction_table_window_context.tsv`). Their classes: Well supported 0, Limited 3, Ambiguous (repeat) 1, Weak 254, Not supported 269. Of the 269 not supported, 6 had at least one fragment by the standard count; the others arose from single reads, contigs or soft-clip consensus sequences that, once each read competes for its best placement, are explained at least as well by an unrearranged repeat copy.

| Sample | Events with ≥ 1 fragment before rescue | after rescue | after rescue, ≥ 2 fragments |
| :--- | :---: | :---: | :---: |
| S01 | 24 | 27 | 0 |
| S02 | 18 | 21 | 0 |
| S03 | 22 | 26 | 1 |
| S04 | 17 | 22 | 1 |
| S05 | 21 | 17 | 2 |
| S06 | 24 | 29 | 0 |
| S07 | 30 | 36 | 0 |
| S08 | 16 | 20 | 0 |
| S09 | 32 | 40 | 0 |
| S10 | 24 | 29 | 0 |

The earlier report's candidates with at least one MAPQ ≥ 20 split read in S01–S10 (179 events; the basis of the per-sample numbers in our 2026-10-14 message) resolve as follows: 175 Weak, 1 Limited (EV0623 = J0786), 2 Not supported, 1 without a counterpart in the rescued evidence.

**Caller concordance** for the 527 target-region junctions:

| Callers with a two-end match | Junctions |
| :--- | :---: |
| DELLY + GRIDSS + SvABA | 1 (J0786, S05) |
| Manta + GRIDSS | 1 |
| GRIDSS only | 28 |
| SvABA only | 1 |
| Manta and DELLY together | 0 |
| none | 496 |

For your question on Manta and DELLY: among the earlier Section 7.3 events, the S05 chr21 junction (EV0623) is called by DELLY only (unfiltered per-sample output) and no event is called by both. Agreement between callers run on the same reads is algorithmic agreement, not independent validation, and the absence of a call is not evidence of absence: two of the four junctions in Section 7.4, the S04 telomere junction and the S05 fold-back inversion, are called by no caller; the telomere junction could not be placed by any caller because its partner sequence occurs at every chromosome end. 11 of the 258 events with at least one fragment in S01–S10 also have a fragment in a Mutu I or DG75 sample; none of the four events in Section 7.4 does.

![Figure 5 — before vs after rescue](figures/Fig5_before_after_and_classes.png)

*Figure 5. Left: number of target-region junctions with support in each sample before rescue (standard count, MAPQ ≥ 20) and after rescue (≥ 1 and ≥ 2 fragments); a junction present in several samples is counted in each, and the bars show net numbers, not which junctions were gained or lost (per junction: `tables/junction_by_sample_long.tsv`). Right: classes of the 527 target-region junctions on a log scale, with the junction IDs of the rare classes.*

### 7.6 The depth step at base-pair scale

| Sample | Step position (start of 200 bp bin) | 95% CI | Depth ratio left / right of the step |
| :--- | :---: | :---: | :---: |
| S04 | 114,605,200 | 114,605,000–114,605,400 | 0.98 / 0.43 |
| S05 | 114,605,800 | 114,605,600–114,606,000 | 0.83 / 0.37 |

The ratios are relative to the S01–S03 mean after normalisation to each sample's own proximal flank (chr11:114.40–114.58 Mb), so they describe the step itself; S01–S03 and GM12878 show no step at these positions. The S04 telomere junction (114,605,216–218) and the S05 fold-back inversion (114,605,971–973) each lie in their own sample's interval; the two steps are about 600 bp apart and their CIs are adjacent without overlapping.

### 7.7 chr11 B-allele fraction

| Segment | Sites | S01–S03: fraction homozygous | S04 | S05 | S04 median abs(BAF − 0.5) | S05 median abs(BAF − 0.5) |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: |
| chr11p 0–50 Mb | 34,495 | 0.000 | 0.005 | 0.002 | 0.100 | 0.093 |
| proximal 99.0–114.5 Mb | 11,363 | 0.000 | 0.004 | 0.003 | 0.100 | 0.159 |
| distal 114.7–135.0 Mb | 16,650 | 0.000 | 0.978 | 0.980 | 0.500 | 0.500 |

"Fraction homozygous" counts sites heterozygous in all of S01–S03 at which the clone shows BAF < 0.1 or > 0.9. Distal of the step both clones retain a single chr11q homolog — the same one in both (8,004 sites homozygous in both clones; the same allele retained at 99.78%). If the two clones had lost different homologs, they would retain opposite alleles at almost every one of these sites, so agreement near 100% identifies the same homolog; the agreement holds along the whole segment (≥ 98.96% in each of eleven 2 Mb windows from the step to the end of 11q), so no part of the distal arm was lost from different homologs in the two clones. The 18 discordant sites (0.22%) are scattered single sites, not a contiguous segment; all are SNVs of the sub-lineage shared by S02, S03 and S05 (Section 8.3), retained in S05 and never present in S04. Which parent this homolog came from is not known. The comparison covers the distal segment only; on which homolog the graded proximal imbalance of S05 lies was not tested. In S05 the proximal imbalance increases toward the step (Figure 6), together with the 1.48× proximal depth reported earlier. GM12878 (S06–S10) shows no homozygous segment on chr11 (all fractions ≤ 0.007).

![Figure 6 — chr11 B-allele fraction](figures/Fig6_chr11_BAF.png)

*Figure 6. BAF at sites heterozygous in the cell line's no-Dox baseline (DLD1: all of S01–S03; GM12878: S06), shown for every reported sample S01–S10. Because sites are selected for allele fraction 0.25–0.75 in the baseline, the baseline panels (S01–S03, S06) are bounded by that range by construction; S04, S05 and the GM12878 clones S07–S10 are free to deviate from it. Vertical line: start of the target region.*

### 7.8 Clone-specific events at the locus

DLD1, target region: 39 junctions are present in S04 or S05 and in none of S01–S03 (22 in S04, 17 in S05); 3 of them have ≥ 2 fragments (the three +Dox junctions of Section 7.4) and 36 have 1 fragment; 10 of the 39 have at least one underpowered baseline sample. No junction is present in both +Dox clones. GM12878, target region: 122 junctions are present in one clone and absent from S06 (S07 36, S08 20, S09 40, S10 29); each has a single fragment. Per-sample labels are in `tables/junction_table_target.tsv` (columns `DLD1_dox_comparison`, `GM12878_clone_specific`, `S0n_state`).

In the ± 50 kb flanks, 13 further junctions have two or more fragments (Figure 3; `tables/junction_table_window_context.tsv`). Two of them are clone-related: J0902 (S05 only, 2 fragments; chr11:114,629,853–858 joined to an L1PA3 LINE element on chrX, a repeat with hundreds of equally good placements, hence "ambiguous") and J0591 (S04 only, 2 fragments; chr11:114,602,525–528 joined to an AluSx3 SINE element on chr16, inside the low-complexity segment that no sample covers with MAPQ ≥ 20 reads, so no baseline sample can be called absent). Neither lies inside a depth-step interval. The two well-supported flank junctions, J1365 (chr11:114,590,721 to chr18:64,754,327) and J0533 (chr11:114,590,798–799 to chr14:83,887,907–908; called by Manta, DELLY and GRIDSS), are present in all five GM12878 samples with 3–10 fragments each and absent from all DLD1 samples; being shared by the population and every clone, they are a feature of the GM12878 line rather than of any clone. The other flank events with two or more fragments are shared across cell lines or occur in the background panel and carry no clone-specific signal.

### 7.9 Genome-wide DLD1 clone comparison

Every DELLY site found in any of S01–S05 (12,953 DELLY events after merging) was re-genotyped in all five samples and labelled with the criteria of Section 5 (`genome_wide_DLD1_comparison/dld1_genome_wide_events.tsv`):

| Label (column `comparison_S04_S05_vs_S01_S03`) | Events |
| :--- | :---: |
| `dox_specific` | 13 |
| `dox_specific_underpowered` | 9 |
| `shared_with_baseline` | 12,323 |
| `baseline_only` | 436 |
| `not_present` | 172 |

Most events are shared with the baseline because they are germline structural variants of the DLD1 line. The 13 +Dox-specific events (`genome_wide_DLD1_comparison/dld1_dox_specific_events.tsv`):

| Event | Type | Position (GRCh38) | Size | Clone | Variant / reference reads | Allele fraction | DELLY somatic filter | Also called by Manta | Gene at a breakend |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |
| DEL00003167 | rearrangement junction (deletion-type orientation); no copy loss of the 888 kb interval; ~3.4 kb gain at chr4:157,172,208–~157,175,500 | chr4:156,284,227–157,172,208 | 888 kb between breakends | S05 | 22 / 16 | 0.58 | PASS | yes | — |
| DUP00005675 | duplication (read depth 1.46× and allele ratio 1:2 in S04) | chr8:46,747,734–47,074,827 | 327 kb | S04 | 8 / 27 | 0.23 | not passed | yes | — |
| INV00010675 | inversion | chr18:22,608,248–22,663,024 | 55 kb | S04 | 15 / 23 | 0.39 | not passed | yes | — |
| INV00004547 | inversion | chr6:115,967,399–115,967,894 | 495 bp | S04 | 14 / 34 | 0.29 | PASS | no | FRK |
| INV00003974 | inversion | chr5:151,020,586–151,020,672 | 86 bp | S04 | 11 / 14 | 0.44 | not passed | no | GPX3 (one end) |
| DEL00002181 | deletion | chr3:112,468,874–112,468,953 | 79 bp | S04 | 14 / 14 | 0.50 | not passed | yes | BTLA |
| DEL00010676 | deletion | chr18:23,100,723–23,100,779 | 56 bp | S05 | 13 / 14 | 0.48 | not passed | no | — |
| INV00005577 | inversion | chr8:23,068,166–23,068,210 | 44 bp | S05 | 4 / 21 | 0.16 | not passed | no | TNFRSF10B |
| DEL00000595 | deletion | chr1:188,196,512–188,196,555 | 43 bp | S04 | 22 / 14 | 0.61 | PASS | yes | — |
| INV00012183 | inversion | chr22:26,292,974–26,293,013 | 39 bp | S05 | 5 / 47 | 0.10 | not passed | no | SEZ6L |
| DEL00000891 | deletion | chr2:11,618,103–11,618,135 | 32 bp | S04 | 31 / 15 | 0.67 | not passed | no | GREB1 |
| DEL00004257 | deletion | chr6:33,234,198–33,234,229 | 31 bp | S05 | 15 / 9 | 0.63 | PASS | yes | — |
| DEL00004287 | deletion | chr6:40,563,333–40,563,359 | 26 bp | S04 | 19 / 7 | 0.73 | PASS | yes | LRFN2 |

Read counts are DELLY's (variant pairs plus junction reads; Section 5). In every listed event each of S01–S03 has 0 variant reads and at least 11 reference reads. Allele fractions near 0.5 are what a heterozygous event in a clonal population gives; the two inversions at 0.10 and 0.16 rest on 4–5 variant reads and are the weakest entries. The 9 underpowered events (7 in S04, 2 in S05) have no variant read in S01–S03, but at least one baseline sample has fewer than 10 reference reads. DELLY's own somatic filter passes 7 events in total; 5 of them are among the 13. None of the 13 events is present in both +Dox clones, and none lies on chr11: the chr11q loss and its junctions (Sections 7.4–7.7) are not represented in the caller output, which is why the locus was analysed with read rescue. The two largest events were checked against read depth (MAPQ ≥ 20, each sample normalised to its own flanks and then to the mean of S01–S03) and against the allele fractions of SNPs heterozygous in S01–S03. The chr8 duplication is confirmed: depth in S04 is 1.46× across the whole interval, with the step at both breakends, and heterozygous SNPs inside it move from 1:1 to 1:2, as a single-copy gain in a clonal population gives. The chr4 event is not a deletion: depth over the 888 kb interval in S05 is 1.04× and all 358 heterozygous SNPs inside it remain heterozygous. Its junction is real (8 independent fragments in S05, none in S01–S04), and depth rises about 1.5× over the ~3.4 kb immediately beyond chr4:157,172,208; the other end of this segment lies in a LINE-1 (long interspersed nuclear element 1) repeat and could not be placed with short reads. The type column of the delivered table keeps DELLY's original call; the depth, allele-balance and junction-read tables of this check are in `supporting_analyses/SV_verification/`. The other 11 events were not read-rescued or curated individually; they are candidates for inspection in the whole-genome BAM files and for PCR.

### 7.10 Genome-wide SNV and indel calls

Joint calling of S01–S10 produced 10,104,360 records after splitting multi-allelic sites (9,477,203 with QUAL ≥ 30), in `genome_wide_SNV_indel_VCF/`. The chr11 records of this call set are identical in number (451,844 at QUAL ≥ 30) to the chr11 call used for the BAF analysis, which was made with the same parameters. Per-sample genotypes, allele depths (AD: reads supporting each allele) and read depth (DP: number of reads at the site after the MAPQ and base-quality thresholds, not corrected for overlapping mates) are in the FORMAT fields; a summary from `bcftools stats` is in `S01_S10.genome.QUAL30.stats.txt`. These calls were not filtered further and were not used to call clone-specific small variants; they are delivered as the requested SNV/indel files.

### 7.11 PCR primers for confirmation

For two of the 17 junctions with two or more fragments, a primer pair was found in which one primer crosses the junction (it cannot prime unrearranged DNA) and the other binds exactly one reference site, at the expected flank (`tables/primer_design.tsv`, `tables/primer_binding_sites.tsv`):

| Junction | Left primer | Right primer | Product on the junction template | Note |
| :--- | :---: | :---: | :---: | :---: |
| J0786 (S05, chr21) | `TGCATCGTGGGTGACAATACTAA` (crosses the junction) | `GCTGCATTGCCAATGAACAGTA` | 165 bp | A second pair, `GTAATCACGGTTACCCAGTGGAT` / `AGCAATAGTTCTCAACAGTGGGT` (274 bp), has each primer at one reference site, on the two sides of the junction, so it amplifies only across the junction |
| J0922 (S04, telomere) | `GCCTCACAGGTAACCACTGTAA` | `ACCCTAACCCTAACCCTAACGA` (crosses the junction) | 234 bp | The right primer is mostly telomeric repeat and will anneal at every telomere; only its 3′ end is chr11 sequence, so background is expected |

For the other 15 junctions, including the S05 fold-back, every candidate pair had at least one primer with more than one reference site, or the template is an inverted duplication on which two-primer PCR is not applicable; the reason is given per junction. These are in-silico designs; no PCR was run.

---

## 8. Conclusions, Limitations and Next Steps

### 8.1 Conclusions by sample group

| Sample group | Conclusion | Basis |
| :--- | :--- | :---: |
| S04 (+Dox clone 3) | Terminal loss of one chr11q homolog from 114,605,216–218, with telomeric repeat joined to the retained chr11 end | 4 fragments, 7 read ends crossing the junction; junction inside the depth-step CI; distal loss of heterozygosity at 97.8% of sites (Sections 7.4, 7.6, 7.7) |
| S05 (+Dox clone 4) | Loss of the same homolog distal of 114.606 Mb; a fold-back inversion at the step and a chr21 junction 1.2 kb distal; proximal gain with graded imbalance | 2 + 3 fragments; fold-back inside the depth-step CI; distal loss of heterozygosity 98.0%; proximal BAF imbalance 0.159; chr21 junction called by three callers (Sections 7.4–7.7) |
| S01–S03 | Baseline for the +Dox clones; carries none of the three +Dox junctions | ≥ 11 reference fragments per sample at each junction, 0 supporting (Section 7.4) |
| S06–S10 | No clone-specific junction above one fragment at the locus; a GM12878-wide junction pair 13 kb proximal of the target | Sections 7.5, 7.8 |
| Locus overall | Most single-read junctions in this 22 kb repeat cannot be distinguished from repeat-copy misplacement; 269 of 527 are not supported once reads compete for placement | Sections 7.1, 7.5 |
| DLD1 genome-wide | 13 clone-private +Dox-specific candidate SVs (8 in S04, 5 in S05), none shared by the two clones | Each present (≥ 3 variant reads) in one +Dox clone and absent (0 variant reads, ≥ 11 reference reads) in all of S01–S03; 5 pass DELLY's somatic filter, 7 are also called by Manta (Section 7.9) |

### 8.2 Practical reading

This section explains, for each sample group, what each result is and what it means, in more detail than Section 2; the notes below the table state the interpretation limits of the S04 and S05 rows.

| Sample group | Samples | Result | What it means |
| :--- | :---: | :---: | :---: |
| **DLD1 +Dox clone 3** | S04 | The broken end is resolved at base-pair level: chr11 sequence ends at 114,605,216–114,605,218 and continues directly into the telomeric repeat (TTAGGG)n. The position lies inside the depth step, and distal of the break one of the two chr11 homologs is lost to the end of 11q. | A terminal deletion of 11q whose broken end carries telomeric repeat, the pattern of a broken chromosome end capped by telomere sequence. A junction-specific PCR primer pair is provided. |
| **DLD1 +Dox clone 4**: overall | S05 | Distal of the step, the same homolog is lost as in S04 (Section 7.7); the proximal side shows a graded allelic imbalance. Two junctions lie near the step, both absent from S01–S03 (next two rows). | A more complex end than in S04; the combination is consistent with repeated breakage and rejoining at this site. |
| **DLD1 +Dox clone 4**: fold-back inversion | S05 | A 376 bp fold-back inversion inside the depth step: chr11:114,605,973 joined head-to-head to 114,606,349. | The 376 bp is not inserted foreign sequence: chr11 is read up to 114,605,973 and then continues backwards from 114,606,349, i.e. the broken end has folded back on itself. That homolog therefore carries no chr11 sequence beyond this point, which matches the distal loss. |
| **DLD1 +Dox clone 4**: chr21 junction | S05 | A junction from chr11:114,607,192 to chr21:20,648,241, 1.2 kb distal of the step; the event called by three of the four callers. | Only the chr11 sequence up to 114,607,192 is joined to chr21; the remaining ~20 Mb of distal 11q is not part of it and is still present on only one homolog, so the junction does not contradict the loss. How it relates to the fold-back is not resolved; Section 7.4 sets out three possible readings. |
| **DLD1 no-Dox population and controls** | S01–S03 | None of the three S04/S05 junctions is present (each with at least 11 reference fragments examined); S03 carries one 2-fragment local event of its own. | They remain a valid baseline: the S04/S05 junctions are specific to the +Dox clones. |
| **GM12878 population and clones** | S06–S10 | No clone-specific junction above one fragment at the locus. All five GM12878 samples share two well-supported junctions 13 kb proximal of the target (chr11:114,590,721 and 114,590,799), absent from every DLD1 sample. | No clone-specific rearrangement at the locus above the detection limit; the shared junction pair is a feature of the GM12878 line itself, not of any clone. |
| **All 527 target-region junction candidates** | S01–S10 | After read rescue, 4 have two or more independent fragments; 254 rest on a single fragment; 269 are explained at least as well by an unrearranged copy of the repeat. | The locus is a 22 kb tandem repeat; most single-read junctions cannot be distinguished from repeat-copy misplacement. |
| **DLD1, genome-wide** | S01–S05 | 13 structural variants elsewhere in the genome are present in S04 or S05 and absent from all of S01–S03 (8 in S04, 5 in S05); none is present in both +Dox clones and none lies on chr11. 7 of the 13 are also called by Manta. | Candidate clone-specific events outside the locus, each private to one clone; none points to a second rearrangement shared by the two +Dox clones. |

**Notes on the table**

1. **Row 1.** "Telomeric repeat" is a sequence description: the reads show only (TTAGGG)n beyond the junction, a sequence present at every chromosome end, so short reads cannot say whether new telomere was added at the break or the break joined an existing chromosome end (Section 8.3). The orientation of the repeat is the one found on a natural chromosome end.
2. **Row 2.** "Consistent with repeated breakage and rejoining" is our interpretation of three measurements together (a fold-back inversion at the step, the proximal gain and graded imbalance, the distal loss); it is not a measurement in itself.
3. **Rows 1–4.** Each event rests on 2–4 independent DNA fragments in one sample. The highest claim level any junction reaches in this report is "limited, coherent with the copy-number step" (Section 5): the sequence evidence is exact, but no independent assay has been run. PCR across each junction is the direct confirmation (Section 8.4).

### 8.3 Limitations

* **Partner of the S04 junction.** The reads contain only (TTAGGG)n beyond the junction. They cannot distinguish new telomere added at the break from a join to the telomere of another chromosome end. In every supporting read the sequence beyond the junction is telomeric repeat up to the end of the read or of the DNA fragment.
* **Support level.** The three +Dox junctions rest on 2–4 independent fragments each in one sample (about 10–15 fragments cover each position). None reached the "Well supported" class, and re-aligning the same sequencing data, however carefully, is not an independent biological test; the coherence tests of Section 5 (depth step and allele fractions) comes from the same sequencing run. PCR across each junction is the direct confirmation (Section 8.4).
* **Detection limit.** With 10–28 fragments spanning a position, a junction carried by a fraction 0.25 of the spanning DNA fragments yields ≥ 3 supporting fragments with probability about 0.5–0.97, and one at 0.10 with probability about 0.1–0.5 (per-sample values: `tables/junction_by_sample_long.tsv`). Junctions in a minority of cells may be missed; "absent" means "not found above roughly this level".
* **Repeat array.** Within the 42 bp-unit array, a rearrangement that joins two near-identical copies (an expansion or contraction of the array) produces a sequence identical to the reference and cannot be seen at the junction level; such events would show only as a depth change.
* **Window.** Junctions are examined when at least one end lies within the target ± 50 kb; events with both ends farther away are covered only by the genome-wide comparison of Section 7.9, which relies on caller output without read rescue.
* **Genome-wide comparison.** Section 7.9 rests on DELLY re-genotyping and per-sample Manta calls; its read counts are DELLY's and are not de-duplicated fragments, and its events, apart from the two largest, were not curated individually. SNV/indel calls were made for all chromosomes but not compared between clones.
* **Sub-lineages within DLD1.** S04 belongs to a different sub-lineage of the DLD1 population from S02, S03 and S05. Among single-nucleotide variants (SNVs) absent from all five GM12878 samples, 84,480 on the autosomes are heterozygous in S02, S03 and S05 and absent from S04, whereas S04 carries 121,783 that none of the other three clones carries; in the population S01 these two sets are present at median allele fractions of about 0.10 and 0.05, i.e. in roughly one fifth and one tenth of its cells. (This comparison was made only to establish how the clones are related; no clone-specific SNV calls are reported. Counts per chromosome and the 18 discordant sites of Section 7.7: `supporting_analyses/DLD1_sub_lineage/`.) The baseline therefore contains no clone from the S04 sub-lineage, only the population, in which a variant of that sub-lineage is expected at an allele fraction of about 0.05, usually 0–1 supporting reads at the depth of S01. Events labelled +Dox-specific in S04, including the chr8 duplication of Section 7.9, can therefore not be distinguished from variants already present in the S04 sub-lineage before doxycycline treatment. S05 shares its sub-lineage with the control clones S02 and S03 and is not affected in this way. The two chr11q losses occurred in different sub-lineages, which supports that they are independent events (Section 7.7).
* **Clone comparison.** The +Dox clones are n = 2, and there is no doxycycline-only control without induction of Epstein–Barr virus nuclear antigen 1 (EBNA1) (as stated in the earlier report).

### 8.4 Suggested confirmation and next steps

| Open question | Suggested experiment |
| :--- | :--- |
| Is the S04 junction real, and is it a telomere added at the break? | PCR with the J0922 pair (Section 7.11) or with a chr11 primer proximal of 114,605,216 and a telomere-repeat primer (single-telomere-length-type assay) on S04, S01–S03 and S05 DNA; FISH with a probe proximal of 114.6 Mb and a telomere probe |
| Are the S05 fold-back inversion and chr21 junction real? | Junction-spanning PCR on S05 and S01–S04 DNA with the J0786 pairs (Section 7.11); the fold-back needs a design outside the duplicated segment (`sequences/junction_and_reference_templates.fa`) |
| Does S05 show further fold-back or rejoining products along proximal 11q? | Copy-number analysis of 11q at higher resolution (for example droplet digital PCR at several proximal positions); deeper short-read sequencing of S05 would raise junction detection along the proximal gain |
| Is the same chr11q homolog targeted in further +Dox clones? With two clones, the same homolog is expected by chance in half of cases; a preference would require more independent clones (for example, all of six clones on the same homolog has a probability of about 3% by chance). | BAF and loss-of-heterozygosity analysis of additional independent +Dox clones (existing pipeline; no new method needed) |
| Are the genome-wide +Dox-specific candidates real? | Inspection of the listed events in the whole-genome BAM files in IGV, then PCR for those of interest (Section 7.9) |

---

## 9. Deliverable Files

All paths are relative to this delivery folder. This report is the entry point.

| File or folder | Content |
| :--- | :--- |
| `Chr11_Locus_and_DLD1_Clone_Analysis_Report_<MMDD>.md` | This report |
| `tables/junction_table_target.tsv` | One row per junction with an end in the target region: breakpoints and uncertainty intervals, type, microhomology and inserted sequence, evidence class, claim level, ambiguity reason, alternative partner placements, flank uniqueness, repeat annotation, columns for Manta, DELLY, GRIDSS and SvABA, earlier event ID, and per sample S01–S10 the before/after fragment counts, reference fragments, allele fraction with 95% CI and present/absent/underpowered state |
| `tables/junction_table_window_context.tsv` | Same columns for junctions with ends only in the ± 50 kb flanks |
| `tables/junction_by_sample_long.tsv` | Every junction × all 18 samples (S11–S18 marked as background panel): counts, allele fraction, detection probabilities |
| `tables/caller_concordance_by_sample.tsv` | Every junction × sample × call set (Manta, DELLY, GRIDSS, SvABA): match class, smallest two-end distance, record IDs |
| `tables/structural_coherence_tests.tsv` | For every junction with ≥ 2 fragments in a sample: the depth-step interval, the breakpoint intervals inside it, and the loss-of-heterozygosity fraction (the two tests of Section 5); column `coherent` is 1 when `criterion1_in_step_ci` and `criterion2_LOH_ge_0.90` are both 1 |
| `tables/key_junction_read_evidence.tsv` | For every junction with ≥ 2 fragments: every supporting read and all of its original alignment records (flag, position, CIGAR string — the per-base alignment description, MAPQ, SA/XA tags, duplicate flag, alignment scores, strand, query interval, read sequence, fragment ID, and the segment overlap or gap of the read end) |
| `tables/rescue_source_summary.tsv` | Rescued read pairs per sample by source; completeness |
| `tables/evidence_filter_impact.tsv` | Count of every evidence item removed and the reason |
| `tables/primer_design.tsv`, `tables/primer_candidates.tsv`, `tables/primer_binding_sites.tsv` | Selected junction-spanning primer pair per junction with ≥ 2 fragments (or the reason none is specific), every candidate pair with its reference match counts, and the location of every reference site of the selected primers |
| `tables/client_request_coverage.tsv` | Every question and deliverable in your request: status, report section, delivered file(s), limits |
| `tables/external_input_status.tsv` | Inputs not available to this analysis and their effect on interpretation |
| `tables/data_inventory.tsv` | Every delivered file with its size |
| `methods/header_sanitization_audit.tsv`, `methods/delly_conversion_audit.tsv` | Removal of local paths from BAM/VCF headers, and DELLY BCF-to-VCF conversion, each with record counts and checksums before and after |
| `figures/` | Figures 1, 2, 3, 5 and 6 (`Fig1_locus_overview`, `Fig2_rescue_sources`, `Fig3_key_events_summary`, `Fig5_before_after_and_classes`, `Fig6_chr11_BAF`), PNG and PDF; `figure_manifest.json` lists every figure file |
| `figures/junction_views/` | Figure 5: one base-level view per junction with ≥ 2 fragments (17 junctions; PNG and PDF) |
| `locus_BAM/` | For S01–S10: every alignment record of every rescued read pair (`<sample>.locus_rescue.bam` and index), viewable in the Integrative Genomics Viewer (IGV); tag `ZR` gives the rescue source |
| `sequences/junction_and_reference_templates.fa` | Junction and reference sequences (400 bp each side) for every junction; basis for primer design |
| `sequences/local_assembly_contigs.fa` | All locally assembled contigs |
| `supporting_analyses/depth_boundary/` | `depth_bins_200bp.tsv`, `boundary_estimates.tsv` |
| `supporting_analyses/chr11_BAF/` | `baf_segments.tsv`, `retained_homolog_S04_vs_S05.tsv`, `baf_sites.tsv.gz` |
| `supporting_analyses/SV_verification/` | Check of the two largest +Dox-specific events of Section 7.9: read depth in bins (`*_depth_bins.tsv`, `*_depth_summary.tsv`, `*_fine_depth_*.tsv`), allele balance at baseline-heterozygous SNPs (`*_allele_balance.tsv`), and reads at the chr4 junction (`DEL00003167_junction_reads.tsv`, `DEL00003167_S05_split_records.tsv`) |
| `supporting_analyses/DLD1_sub_lineage/` | SNVs shared between DLD1 clones (`clone_snv_sharing_autosomes.tsv`, `clone_snv_sharing_by_chrom.tsv`) and the 18 distal-11q sites at which S04 and S05 retain different alleles (`distal11q_discordant_sites.tsv`) (Section 8.3) |
| `genome_wide_SV_VCF/` | Per sample S01–S10: Manta `diploidSV` unfiltered and PASS, Manta `candidateSV`, DELLY unfiltered and PASS (with index); `vcf_inventory.tsv`. `locus_GRIDSS_SvABA/`: the GRIDSS (18 samples) and SvABA (S01–S05, filtered and unfiltered) calls for the analysis window |
| `genome_wide_SNV_indel_VCF/` | `S01_S10.genome.raw.vcf.gz` (all calls), `S01_S10.genome.QUAL30.vcf.gz` (QUAL ≥ 30), indexes, `S01_S10.genome.QUAL30.stats.txt` |
| `genome_wide_DLD1_comparison/` | `dld1_genome_wide_events.tsv` (every re-genotyped DELLY event in S01–S05 with per-sample counts, states, comparison label, Manta column and genes), `dld1_dox_specific_events.tsv` (the +Dox-specific subset), `summary_counts.tsv`, `dld1_S01_S05.delly_genotyped.vcf.gz`, `somatic_samples.tsv` |
| `md5sums_delivery.txt` | Checksums of every file in this folder |
| Whole-genome BAM files | Delivered earlier; not repeated in this package |

---

**Zhen Gao, PhD**\
Principal Bioinformatics Scientist, Athenomics
