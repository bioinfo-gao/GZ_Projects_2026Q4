# v1_20261027: Chinese reading copy of the English report as of c32b1d8 (before the final 1027 review edits)

- What: `P45_Chr11_Locus_and_DLD1_Clone_Analysis_Report_1027_CN.md` as it was at 17:17 on 2026-10-27, a translation of the
  English report at commit c32b1d8.
- Why superseded: the user reviewed the English 1027 report (审阅完毕, 2026-10-27) after these edits: header/signature hard
  breaks; Section 2 and Section 8.2 bodies swapped with a strengthened list; S05 split into three rows; Methods UCSC table
  names; Section 7.4 three readings of the chr21 junction. The user asked for the CN copy to be regenerated and to replace
  the old file under the same name.
- How: `update_CN_to_final_1027.py` (snapshot here) applies the same 7 changes to this file. The rest of the translation
  is unchanged, because it already mirrored the unchanged English text line by line.
- Compare: `diff P45_Chr11_Locus_and_DLD1_Clone_Analysis_Report_1027_CN.md ../../../docs/P45_Chr11_Locus_and_DLD1_Clone_Analysis_Report_1027_CN.md`
