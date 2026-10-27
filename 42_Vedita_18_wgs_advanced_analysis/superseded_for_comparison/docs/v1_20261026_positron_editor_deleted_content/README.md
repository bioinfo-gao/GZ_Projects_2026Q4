# v1_20261026 — two internal CN docs damaged by a Positron WYSIWYG Markdown editor (found 2026-10-26)

Both files were rewritten by the IDE (Positron local history source "Workspace Edit"), and the damaged versions were then
committed unnoticed. After normalising whitespace and table padding, the deletions below are the ONLY content changes.

| file | damaged at | lost | committed in |
| :--- | :---: | :---: | :---: |
| docs/P42-updated1/P42_analysis_rationale_and_process_20261024_CN.md | 2026-10-24 19:28:44 | §7 "已建立" line cut after "…crosswalk；C"; "未建立" and "合理的后续证据" paragraphs deleted; trailing two-space hard breaks of the header stripped | 5cd4863 |
| docs/P42-updated1/P42_Client_Request_Completion_1024_CN.md | 2026-10-25 10:21:02 | key point 5 (source-data inventory, 907.44 GB BAM …) deleted; header/signature hard breaks stripped | 9291bee |

- Suspect: "墨览 Markdown" (fengshihao.molan-markdown), default custom editor for *.md since 2026-10-24 07:31; uninstalled
  2026-10-26. Office Viewer also registers a *.md editor; attribution not provable.
- Restored 2026-10-26 from: rationale = Positron local-history snapshot saved 2026-10-24 18:56 (last save before the
  damage; copy kept here as `*.pre_damage_positron_history.md`); completion = git c39555b.
- Files: `*.damaged.md` = state on disk (and in git) before the restore.
- Compare: `diff P42_analysis_rationale_and_process_20261024_CN.damaged.md ../../../docs/P42-updated1/P42_analysis_rationale_and_process_20261024_CN.md`
