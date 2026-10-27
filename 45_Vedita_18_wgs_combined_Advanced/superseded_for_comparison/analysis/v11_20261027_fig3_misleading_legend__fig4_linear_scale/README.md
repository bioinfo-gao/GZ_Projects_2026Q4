# v11_20261027 — Fig3 legend misleading / zero categories; Fig4 rare classes invisible, no cell lines

- Produced: 2026-10-25 16:34. Original delivered files, copied before replacement.
- Fig3 defects: green legend "in W + mate/SA/XA outside" was wrong (13,040 of 13,093 such pairs in S01 have the mate IN
  the window); three legend entries (b, c, e) had 0 pairs and no visible bar; y labels sample number only; no note that
  counts are not depth-normalised.
- Fig4 defects: right panel on a linear scale made the 3 Limited and 1 Ambiguous junctions (the important ones) nearly
  invisible; x labels sample number only; "0" labels sat on the x axis (no collision check was run on Fig4 before).
- Fix (v12, user approved 2026-10-27): Fig3 two segments (a; d only) with an assertion per sample that a + d-only = total
  and b/c/e-only = 0, footnote for b/c/e, cell-line labels, depth note; Fig4 cell-line tick labels, right panel log scale
  with junction IDs for classes with 1-5 events, "0" labels offset, collision check on both figures. Captions rewritten
  in all three report copies. No numbers changed.
- Contents: figures/ (old PNG/PDF), scripts/11_figures.py (snapshot before the change), manifest, md5, old captions.
