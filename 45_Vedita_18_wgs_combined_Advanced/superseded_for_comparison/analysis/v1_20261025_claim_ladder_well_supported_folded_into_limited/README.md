# v1_20261025 — claim ladder folded "Well supported" into "limited junction candidate"

- Produced: 2026-10-25 14:39–14:44 by p45_chain (stage 9x, 10, 11; 10c table copied here as well).
- Defect: plan section 4.2 defines no claim level for the P44 evidence class Well_supported; 9x v1 mapped it to
  limited_junction_candidate, understating two flank junctions (J1365, J0533; 8-10 fragments in all five GM12878 samples).
- Fix (v2): new level well_supported_junction_candidate for Well_supported events that are not structurally coherent.
  No other rule changed. Stage 14/14b outputs did not depend on claim_level and were not re-run.
- Compare: diff scripts/9x_structural_coherence.py ../../../scripts/20261025/9x_structural_coherence.py ;
  diff stage9x_coherence/claim_ladder.tsv ../../../analysis/20261025_p45/stage9x_coherence/claim_ladder.tsv
- These are the original files (moved, not regenerated).
