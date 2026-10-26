# v2_20261025 — SV records' END ignored (pysam does not expose END via rec.info)

- Produced: 2026-10-25 15:50 (stage 10/11 after the claim-level fix) and 16:10 (G2b). Original files, moved.
- Defect: svmatch.pairs() (copied from P44 8_caller_concordance.py) used info.get("END", pos). pysam never puts END in
  rec.info, so every DEL/DUP/INV record became a zero-length event at POS:
  * stage8x: intrachromosomal Manta/DELLY records could not reach a two-end match when the junction's second end is
    near the record END -> concordance for such junctions under-counted (same defect exists in P44 Stage 8);
  * G2b: pos_B = pos_A and size 0 for all DEL/DUP/INV events; Manta matching for those degenerated to one-end tests.
- Fix (v3): end = rec.stop when END is not in rec.info (rec.stop = END for symbolic SV alleles; = POS for BND).
- stage9x did not depend on callers; a copy is kept here only for completeness (not re-run).
- Compare: diff scripts/svmatch.py ../../../scripts/20261025/svmatch.py
