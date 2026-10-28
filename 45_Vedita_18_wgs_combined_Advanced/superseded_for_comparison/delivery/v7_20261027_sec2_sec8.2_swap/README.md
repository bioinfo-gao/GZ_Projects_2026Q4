# v7_20261027: Section 2 and Section 8.2 bodies swapped (this report only)

- Why: the user judged Section 8.2 "Practical reading" (numbered conclusions with evidence) to be what the client most
  wants to see, and asked to swap the bodies of Section 2 and Section 8.2 and to strengthen the list. This applies to
  this report only. The report skill's general structure (Section 2 brief, Section 8.2 detailed) is NOT changed.
- What changed (script `swap_sec2_sec8_2.py`, snapshot in this folder; run on docs/_client_report_source.md):
  - Section 2 now holds the 6-point list. Each point is a bold conclusion, then "What it means", then the original
    Evidence list. Point 6 (genome-wide) got a real conclusion; point 2 (S05) was brought in line with the fold-back and
    chr21 explanations and points to Section 7.4.
  - Section 8.2 now holds the former Section 2 table (8 rows) and its 3 notes. "Rows 1–2" in note 3 became "Rows 1–4"
    because S05 now spans three rows.
  - Both intro sentences were rewritten.
- Checks: every former Section 2 table row and note appears exactly once in the new file; all Section references resolve;
  all 47 file names cited in the report exist in the delivery folder; md5sum -c passes.
- Compare: `diff Chr11_Locus_and_DLD1_Clone_Analysis_Report_1027.before_swap.md ../../../custom_research_report_20261027/Chr11_Locus_and_DLD1_Clone_Analysis_Report_1027.md`
