# P42 chr11 junction read rescue — research plan and review

- **Created:** 2026-10-22
- **Status:** formal analysis authorised by the user; outputs are new P42 products and do not alter the frozen P41 delivery.
- **Scope:** S01–S10 only. This includes S02, S03 and S10. Mutu I and DG75 samples are excluded.
- **Target:** `chr11:114,604,000-114,625,500` (GRCh38/hg38 coordinate system).
- **Canonical sample table:** `scripts/0_sample_table.csv`.

## Objective

Perform the contracted Analysis 1 read-rescue extension of P41: re-extract the fragment sequences supporting candidate junctions in the target, independently remap them against the same GRCh38 reference with all alignments retained, and report whether each candidate's two breakpoint sides are reconstructed by the rescued fragment. This is not manual base-pair curation or local assembly; those were a separately scoped service in the correspondence.

## Inputs and provenance

The authoritative BAM input is `/home/gao/Dropbox/Vedita_Bam/`, which contains mark-duplicate BAMs and indexes copied from P41. On 2026-10-22, all 36 files pass `md5sum -c MD5.txt` when run in that directory; the ten in-scope BAMs pass `samtools quickcheck`.

The starting candidates are P41's per-sample `analysis/region/evidence/target/*.junctions.tsv`. The prior e-mail listed 150 candidates for eight samples but did not specify a machine-reproducible eligibility rule. P42 fixes the rule before analysis: retain a P41 target-end junction only when it has at least one non-duplicate split read whose two P41 segments both have mapping quality (MAPQ) >=20. Under that definition the ten-sample candidate list contains 179 rows: S01 20, S02 14, S03 13, S04 17, S05 12, S06 22, S07 24, S08 13, S09 26 and S10 18. The source evidence remains unchanged in P41.

The source BAM total for the ten samples is approximately 482 GiB (directory sizes measured 2026-10-22). Only the target window is queried; no BAM is copied or edited.

## Analysis steps

1. Materialise the 179-row candidate table and validate that all ten samples, including S02/S03/S10, are present.
2. For each P41 support fragment, extract its sequence from `chr11:114,554,000-114,675,500` without changing the BAM.
3. Remap each extracted sequence with `bwa-mem2 mem -a -Y` to the original GATK GRCh38 analysis-set reference. `-a` preserves alternative placements; this prevents a low-complexity/repeat-derived sequence from being presented as uniquely rescued.
4. Test whether a remapped fragment independently reconstructs both recorded breakpoint ends within 20 bp. Classify every candidate as uniquely reconstructed, reconstructed but ambiguous/low-MAPQ, not reconstructed, or not evaluable because its source sequence was not re-extracted.
5. Produce a per-candidate table, a per-sample summary, a client-facing supplementary report, and a supplementary delivery directory. No conclusion of structural-variant confirmation is made from a single read.

## Pre-registered decision criteria

- A **unique rescued fragment** has alignments to both candidate sides within 20 bp and both matching alignments have MAPQ >=20 after independent remapping.
- A candidate with >=2 unique rescued fragments is labelled `multi_fragment_supported`; one is `single_fragment_supported`; one or more reconstructed fragments without the MAPQ criterion are `ambiguous_reconstruction`; otherwise it is `not_reconstructed` or `not_evaluable`.
- The result describes read-level support and mapping ambiguity, not a validated breakpoint sequence. Manual curation/local assembly is outside this authorised Analysis 1 scope.

## Resources and safety

The re-alignment is bounded to a few hundred extracted reads, uses eight threads, and writes only under `analysis/20261022/` plus the new supplementary delivery directory. It is expected to complete within minutes; no long-running pipeline or tmux watchdog is required. P41 and its 2026-10-09 delivery are read-only inputs.

## Change log

- 2026-10-22 — Created after the user authorised the full P42 workflow and expanded scope to S01–S10.
