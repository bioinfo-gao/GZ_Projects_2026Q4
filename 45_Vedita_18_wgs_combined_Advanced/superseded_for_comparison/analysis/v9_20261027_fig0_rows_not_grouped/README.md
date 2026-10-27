# v9_20261027 — Fig0 rows not grouped; blanks ambiguous; no cell-line header; no background column

- Produced: 2026-10-25 16:37 (Fig0 unchanged until 2026-10-27). Original delivered files, copied before replacement.
- Defects (user review 2026-10-26/27, reading guide §8.3 and §10 item 11):
  1. rows ordered "coherent first, then by ID": junctions of the same pattern were scattered;
  2. a blank cell did not distinguish "0 supporting fragments" from "underpowered (< 10 reference fragments)";
  3. header showed only S01-S10, not the cell lines;
  4. presence in background cell lines S11-S18 was not visible;
  5. no position and no target/flank tag in the row labels.
- Fix (v10, approved by user 2026-10-27): rows grouped A-E by a rule on existing table fields (A case_specific*,
  B baseline_only, C not_present & S06 present, D not_present & S06 absent, E shared_with_baseline); within a group target
  first, then by position; two-level DLD1 / GM12878 header with pop/clone/+Dox; underpowered blanks hatched; S11-S18
  background column; labels carry chr11 position and T/F; caption rewritten in all three report copies.
  No numbers, classes or claim levels changed. Figure number unchanged (renumbering deferred by the user).
- Contents: figures/ (old PNG/PDF), scripts/11_figures.py (snapshot before the change), manifests, md5, old captions.
