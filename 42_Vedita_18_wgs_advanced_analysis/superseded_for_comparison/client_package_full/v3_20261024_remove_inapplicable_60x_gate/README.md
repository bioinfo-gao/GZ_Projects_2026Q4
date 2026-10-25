# P42 client-package interpretation correction — 2026-10-24

This internal archive preserves the complete `custom_research_report_20261024` package generated before the depth-interpretation correction. It was prepared for review and was **not sent externally**.

## What required correction

The archived report and selected result tables used a retrospective `60x` local-depth gate when interpreting C1–C4. That gate was not appropriate for this study because it conflated two distinct quantities:

1. the WGS data-yield context (approximately 25–30x for most samples, as supplied by the project owner); and
2. the more restrictive primary, non-duplicate, non-QC-fail, non-supplementary, MAPQ≥20 autosomal alignment denominator used only in the candidate-comparison audit.

The restrictive denominator is not a service-acceptance or service-quality measurement. A comparator zero within that denominator means only that the candidate was not observed under the stated evidence context. It does not establish absence.

## Corrected replacement

The rebuilt `custom_research_report_20261024/` package is created separately from this archive. It retains the validated BAM/VCF/figure payload unchanged, removes the retrospective numeric depth gate, distinguishes WGS data yield from the comparison denominator, and keeps C1–C4 as candidate observations for which clone- or doxycycline-specific status is not established from the current WGS comparison alone.

The historical package, scripts snapshot, and this explanation remain internal only and are not part of the replacement client package.
