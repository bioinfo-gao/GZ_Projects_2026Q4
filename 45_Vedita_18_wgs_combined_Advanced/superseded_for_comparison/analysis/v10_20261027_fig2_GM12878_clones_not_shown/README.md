# v10_20261027 — Fig2 showed GM12878 only as S06 (its own site-selection baseline)

- Produced: 2026-10-25 16:34. Original delivered files, copied before replacement.
- Defect: Fig2 drew S01-S05 and S06 only. Sites are selected as heterozygous in S06 for GM12878, so the S06 panel is
  bounded to 0.25-0.75 by construction and cannot show LOH; the GM12878 clones S07-S10 (data present in baf_sites.tsv.gz)
  were not drawn, although the report text states "GM12878 (S06-S10) shows no homozygous segment on chr11". Panel labels
  gave the sample number only.
- Fix (v11, user approved 2026-10-27): all ten reported samples drawn (S07-S10 added); panel labels name cell line and
  population/clone (+Dox); assertion that every reported sample has a panel; caption rewritten in all three report copies.
- Contents: figures/ (old PNG/PDF), scripts/11_figures.py (snapshot before the change), manifest, md5, old captions.
