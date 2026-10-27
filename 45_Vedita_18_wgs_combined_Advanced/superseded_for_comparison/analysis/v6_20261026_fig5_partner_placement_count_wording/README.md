# v6_20261026 — Fig5 "AMBIGUOUS PARTNER" line miscounts the placements

- Produced: 2026-10-25 16:37 (stage 11, v6 of the Fig5 layout). Original files, copied before replacement
  (delivery copy and `analysis/20261025_p45/stage11_figures/junctions/` copy were byte-identical; kept once).
- Defect: line 2 of every Fig5 with an ambiguous partner printed
  "B-side sequence (...) has {n_alternative_representations} equivalent placements; the coordinate shown is one of them".
  `n_alternative_representations` counts the alternatives *besides* the coordinate shown, so the text under-counts by one
  and contradicts the report body (J0922: chrX shown + chr10, chr2, chr4 = 4 placements; figure said 3).
- Affected (3 of 17 Fig5): J0922 (3 alternatives -> 4 in total), J1635 (8 -> 9), J1474 (2 -> 3). The other 14 have a
  unique partner and their text is unchanged.
- Fix (v7, 2026-10-26): "has {n} alternative placement(s) besides the coordinate shown ({n+1} equivalent in total)".
  Only the text changed; read rows, colours and layout are unchanged.
- Rerun: `P45_FIG_OUT=analysis/20261025_p45/stage11_figures_rerun_20261026 python scripts/20261025/11_figures.py`
  (new output dir; the original stage11 dir was not rewritten wholesale). Only the 3 affected Fig5 PNG/PDF pairs were then
  copied into `stage11_figures/junctions/` and the delivery `figures/junction_views/`; their manifest entries and
  md5 lines were updated.
- Contents: `figures_junction_views/` (6 old files), `scripts/11_figures.py` (snapshot before the fix),
  `figure_manifest.delivery.json`, `figure_manifest.stage11.json`, `md5sums_delivery.txt` (all pre-fix).
- Compare: `diff scripts/11_figures.py ../../../scripts/20261025/11_figures.py`; open old and new PNG side by side,
  e.g. `figures_junction_views/Fig5_J0922.png` vs
  `../../../custom_research_report_20261025/figures/junction_views/Fig5_J0922.png` (line 2 only).
