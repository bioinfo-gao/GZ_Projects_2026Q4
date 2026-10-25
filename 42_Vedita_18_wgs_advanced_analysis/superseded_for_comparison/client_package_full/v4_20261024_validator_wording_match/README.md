# P42 corrected-package build record — 2026-10-24

This internal record preserves the first build of the corrected P42 package. The package content was created successfully, but the self-validator then failed because it looked for a shorter phrase than the equivalent wording actually written in the depth table.

The failure was a validator string-match defect, not a data, figure, BAM, VCF, or scientific-analysis failure. The build is retained unchanged here with its failure audit. The validator was corrected to check the exact approved depth-definition phrase before the package was rebuilt and fully revalidated.

This directory is internal only and was never sent externally.
