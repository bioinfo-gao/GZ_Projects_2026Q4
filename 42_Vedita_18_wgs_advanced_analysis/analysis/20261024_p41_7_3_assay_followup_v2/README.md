# Superseded analytical output: remote-orientation error

This output was generated during the 2026-10-24 P41 Section 7.3 assay follow-up. It is retained for comparison, but it is **not a valid scientific result**.

The initial virtual-template implementation imposed one remote reference orientation from BND ALT-string handling. Direct inspection of the three C2 split-read QNAMEs instead showed a chr11-forward to chr21-forward sequence transition. As a result, this version reported zero exact junction-spanning 24-mers and consequently made no assay selection.

The corrected non-overwriting rerun evaluates both remote orientations at each coordinate representation and selects the template with maximum exact support in the recovered raw split-read sequences. Its output path is `../20261024_p41_7_3_assay_followup_v3_raw_orientation/`.
