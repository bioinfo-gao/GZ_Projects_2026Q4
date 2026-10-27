# v1_20261026 — delivery report damaged by a Positron WYSIWYG Markdown editor

- What: `custom_research_report_20261025/Chr11_Locus_and_DLD1_Clone_Analysis_Report_1025.md` was rewritten twice by the
  IDE (Positron local history source "Workspace Edit"; 2026-10-26 21:16:54 and 22:21:45). No agent was running at either
  time; the user did not make the edit.
- Damage: every table re-padded (formatting only) and the Figure 5 caption line
  ("*Figure 5 (one of 17; all in `figures/junction_views/`). ...*") deleted. After normalising whitespace and table
  padding, the caption is the ONLY content difference from git c519305.
- Suspect: "墨览 Markdown" (fengshihao.molan-markdown), default custom editor for *.md since 2026-10-24 07:31; uninstalled
  2026-10-26 at the user's request. Office Viewer also registers a *.md editor; which extension wrote it is not provable.
- Fix (2026-10-26): file restored byte-for-byte from git c519305; `md5sum -c md5sums_delivery.txt` passes again.
- Files: `*.damaged_2116_positron_history.md` (21:16 snapshot, copied from Positron local history),
  `*.damaged_2221.md` (state on disk before the restore).
- Compare: `diff <(git show c519305:./45_Vedita_18_wgs_combined_Advanced/custom_research_report_20261025/Chr11_Locus_and_DLD1_Clone_Analysis_Report_1025.md) Chr11_Locus_and_DLD1_Clone_Analysis_Report_1025.damaged_2221.md`
