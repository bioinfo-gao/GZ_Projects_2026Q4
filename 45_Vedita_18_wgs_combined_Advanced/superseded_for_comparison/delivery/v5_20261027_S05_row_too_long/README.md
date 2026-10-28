# v5_20261027: Core Conclusions S05 row too long to read

- Problem: after v3/v4 the single S05 row of the Section 2 table had very long "Result" and "What it means" cells
  (fold-back explanation, chr21 explanation and the pointer to Section 7.4 all in one cell). The user found it painful to read.
- Fix (2026-10-27, user instruction): split into three rows, "overall", "fold-back inversion" and "chr21 junction", one event
  per row. All content is kept and only redistributed. The caller-count sentence moved into the chr21 row's Result cell.
- Compare: `diff Chr11_Locus_and_DLD1_Clone_Analysis_Report_1027.before_row_split.md ../../../custom_research_report_20261027/Chr11_Locus_and_DLD1_Clone_Analysis_Report_1027.md`
