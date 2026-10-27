# v8_20261026 — Fig1 middle panel: S02 and S03 drawn in one colour

- Produced: 2026-10-26 23:29 (v8 = AT-rich block grey-hatched, S01 added). Original delivered files, copied before replacement.
- Defect: the two no-Dox clones S02 and S03 shared one colour (#8B4500) and one legend entry ("S02, S03 (no-Dox clones)"),
  so the reader could not tell which line is which sample (user review 2026-10-26: "完全不对").
- Fix (v9): S02 #B35806 "S02 (no-Dox clone 1)", S03 #1B7837 "S03 (no-Dox clone 2)", separate legend entries; nothing else changed.
- Compare: `diff scripts/11_figures.py ../../../scripts/20261025/11_figures.py`; open `figures/Fig1_locus_overview.png` next to
  `../../../custom_research_report_20261025/figures/Fig1_locus_overview.png` (middle panel only).
