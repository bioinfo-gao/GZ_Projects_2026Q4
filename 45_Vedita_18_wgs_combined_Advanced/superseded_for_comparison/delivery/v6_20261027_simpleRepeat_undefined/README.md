# v6_20261027: UCSC table names simpleRepeat / genomicSuperDups not explained at first use

- Problem: Methods (Section 6, Breakpoint annotation) listed "UCSC hg38 RepeatMasker, simpleRepeat and genomicSuperDups tables"
  without saying what the two table names contain. The user asked whether simpleRepeat is a separate program.
- Fix (2026-10-27, user approved): first mention only, "simpleRepeat (tandem repeats detected by Tandem Repeats Finder)" and
  "genomicSuperDups (segmental duplications)". Later mentions (Section 7.1, Figure 1 caption) are unchanged.
- Compare: `diff Chr11_Locus_and_DLD1_Clone_Analysis_Report_1027.before_table_names.md ../../../custom_research_report_20261027/Chr11_Locus_and_DLD1_Clone_Analysis_Report_1027.md`
- Note: the archived `*.before_table_names.md` is the file as it was on disk at 19:51:48. Two Positron "Workspace Edit" saves
  (19:47:32, 19:51:48) had re-padded every table AND removed the trailing `\` hard breaks from the 8 header lines and the
  signature line. After table re-padding is normalised, those 9 backslashes are the only difference from the source at that
  time. Regenerating from source restored them, which is why a raw diff shows ~366 lines while the only intended change is
  line 157.
