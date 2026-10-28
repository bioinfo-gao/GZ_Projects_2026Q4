#!/usr/bin/env python3
# 1027 review (user instruction, this report only): swap the bodies of Section 2 (Core Conclusions) and Section 8.2
# (Practical reading). Section 2 becomes the numbered list, strengthened with a "What it means" sentence per point;
# Section 8.2 receives the per-sample-group table and its notes. Edits docs/_client_report_source.md in place.
import re, sys

SRC = sys.argv[1]
t = open(SRC).read()

s2_start = t.index("## 2. Core Conclusions\n")
s3_start = t.index("## 3. Key Findings\n")
p2_start = t.index("### 8.2 Practical reading\n")
p3_start = t.index("### 8.3 Limitations\n")
old_s2 = t[s2_start:s3_start]
old_p2 = t[p2_start:p3_start]

# Section 2 body = table + notes (everything after the intro paragraph, up to the closing rule)
m = re.search(r"\n\n(\| Sample group \| Samples \|.*?)\n\n---\n", old_s2, re.S)
table_and_notes = m.group(1)
table_and_notes = table_and_notes.replace("3. **Rows 1–2.** Each event", "3. **Rows 1–4.** Each event")

new_s2 = """## 2. Core Conclusions

This section gives each conclusion in one sentence, what it means for you, and where the evidence is. Section 3 lists the measurements by topic; Section 8.1 sets out the basis of each conclusion; Section 8.2 explains each result per sample group in more detail; Section 8.3 states what the data cannot show; Section 8.4 pairs each open question with an experiment.

1. **S04 (+Dox clone 3): the broken end is resolved — chr11 at 114,605,216–218 joined to telomeric repeat, and distal of it one chr11q homolog is lost to the end of 11q.** What it means: this clone carries a terminal deletion of 11q whose broken end is capped by telomeric repeat. A junction-specific PCR primer pair (Section 7.11) can confirm it on S04 DNA, with S01–S03 as negative controls. Evidence: Section 7.4, `tables/junction_table_target.tsv` (row J0922), `figures/junction_views/Fig4_J0922.png`, `tables/key_junction_read_evidence.tsv`, `locus_BAM/S04_DLD1_c3_dox.locus_rescue.bam`.
2. **S05 (+Dox clone 4): the same chr11q homolog is lost, but the end is different and more complex — a 376 bp fold-back inversion at the step and a chr21 junction 1.2 kb distal of it.** What it means: at the fold-back the broken end has folded back on itself (it is not an insertion); the chr21 junction carries only about 1 kb of chr11 sequence and does not contradict the distal loss; how the two events relate is open (Section 7.4). The chr21 junction can be tested with the J0786 primer pairs (Section 7.11); the fold-back needs a different design (Section 8.4). Evidence: Sections 7.4 and 7.7, rows J0359 and J0786, `figures/junction_views/Fig4_J0359.png` and `Fig4_J0786.png`, `supporting_analyses/chr11_BAF/retained_homolog_S04_vs_S05.tsv`, `tables/caller_concordance_by_sample.tsv`.
3. **S01–S03 remain a valid baseline: none of them carries any of the three +Dox junctions.** What it means: they can serve as negative controls in junction PCR. S04 belongs to a sub-lineage that none of the control clones represents, so events labelled +Dox-specific in S04 cannot be distinguished from variants already present in that sub-lineage before doxycycline treatment (Section 8.3); S05 is not affected. Evidence: Section 7.4 table (reference fragment counts), `tables/junction_by_sample_long.tsv`.
4. **GM12878 clones: no clone-specific rearrangement at the locus above one fragment.** What it means: no event at this locus in S07–S10 calls for follow-up; the junction pair 13 kb proximal of the target, present in all five GM12878 samples, is a feature of the GM12878 line rather than of any clone. Evidence: Sections 7.5 and 7.8, `tables/junction_table_target.tsv` column `GM12878_clone_specific`, `tables/junction_table_window_context.tsv` rows J1365 and J0533.
5. **Most junction candidates at this locus are single reads inside a 22 kb repeat and should not be followed up individually.** What it means: of 527 candidates, 4 have two or more independent fragments after read rescue, and three of these are the +Dox junctions of points 1–2; follow-up effort is best spent on those. Evidence: Sections 7.1 and 7.5, Figure 5, `tables/junction_table_target.tsv` columns `classification` and `claim_level`.
6. **Genome-wide: 13 candidate structural variants are private to one +Dox clone (8 in S04, 5 in S05); none is shared by the two clones and none lies on chr11.** What it means: no second rearrangement common to both +Dox clones was found. Of the two largest events, the 327 kb chr8 duplication in S04 is confirmed by read depth and allele balance, and the chr4 event in S05 is a real junction but not a deletion; the other 11 are candidates for inspection in the whole-genome BAM files and for PCR (Section 8.4). Evidence: Section 7.9, `genome_wide_DLD1_comparison/dld1_dox_specific_events.tsv`, `supporting_analyses/SV_verification/`.

---

"""

new_p2 = ("### 8.2 Practical reading\n\n"
          "This section explains, for each sample group, what each result is and what it means, in more detail than "
          "Section 2; the notes below the table state the interpretation limits of the S04 and S05 rows.\n\n"
          + table_and_notes + "\n\n")

t = t[:s2_start] + new_s2 + t[s3_start:p2_start] + new_p2 + t[p3_start:]
open(SRC, "w").write(t)

# reverse assertions (report skill): no "simple repetition" wording, no apology wording
for pat in [r"the same (answers|ground|findings)", r"repeats the", r"restates the"]:
    assert not re.search(pat, new_s2 + new_p2), pat
print("swapped; Section 2 points:", len(re.findall(r"^\d\. \*\*", new_s2, re.M)),
      "; Section 8.2 table rows:", new_p2.count("\n| **"))
