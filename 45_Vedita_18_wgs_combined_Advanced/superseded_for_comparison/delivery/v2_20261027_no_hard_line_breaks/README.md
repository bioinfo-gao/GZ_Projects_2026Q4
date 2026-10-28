# v2_20261027: report header and signature had no Markdown line breaks

- What happened: in `custom_research_report_20261027/Chr11_Locus_and_DLD1_Clone_Analysis_Report_1027.md` (and its source
  `docs/_client_report_source.md`) the 9 header lines (Report Date to Relation to earlier work) and the signature
  (`**Zhen Gao, PhD**` / title) were separated by single newlines only. Under CommonMark (VS Code/Positron preview, GitHub,
  pandoc) these merge into one paragraph. The old Office Viewer extension showed them on separate lines, so the problem was hidden.
- Fix (2026-10-27, while 1027 was still in internal review): a blank line was inserted after each of the 9 lines in the source,
  and the delivery .md was regenerated with the build script's own `sed s/@@REPORT_DATE@@/2026-10-27/`. Blank lines were
  chosen over `\` or two trailing spaces because they render the same in every Markdown flavour (Python-Markdown and
  Markdown.pl ignore `\`) and also read correctly as plain text. Positron "Workspace Edit" rewrites strip trailing spaces.
  Only the report's md5 line in `md5sums_delivery.txt` changed, and `md5sum -c` passes on all 225 files.
- Content: once blank lines are removed, the new report is line-for-line identical to git HEAD (c32b1d8).
- Files: `*.before_hardbreaks.md` is the delivery file as it was on disk before the fix. It also carries the 17:26 Positron
  "Workspace Edit" table re-padding, which is formatting only and is undone by regenerating from source.
  `_client_report_source.before_hardbreaks.md` and `md5sums_delivery.before_hardbreaks.txt` are the other two pre-fix files.
- Compare: `diff Chr11_Locus_and_DLD1_Clone_Analysis_Report_1027.before_hardbreaks.md ../../../custom_research_report_20261027/Chr11_Locus_and_DLD1_Clone_Analysis_Report_1027.md`
- Revision (same day, user feedback): the blank-line version put a paragraph gap between header fields, which the user rejected
  as ugly; it is kept as `*.blank_line_variant.md`. Final form: a trailing `\` (CommonMark/GFM hard line break) on the 8 header
  lines and on `**Zhen Gao, PhD**`, so the lines stay tight. Any reader that renders this report's GFM tables also honours `\`.
  The build script's Report Date check (`scripts/20261025/12_build_delivery.sh`) now accepts the trailing `\`.
