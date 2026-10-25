# v1 — supplementary delivery as generated 2026-10-25 01:40 (round 1, window W = target ± 5 kb)

- **What this is:** a byte-for-byte copy (`cp -a`) of `custom_research_report_supplementary_20261025/` as built from
  `analysis/20261024_read_rescue/` before the round-2 changes, plus the report source used (`_client_report_source_v1.md`).
  It was never sent to the client and had not yet been reviewed by the user.
- **Why superseded (not a scientific error):** user decisions of 2026-10-25 (`docs/P44-learned-from-P42/P44_learned_from_P42_1025.md` §6–§7):
  1. W extended from ± 5 kb to ± 50 kb (D3) — new data universe, analysis re-run to `analysis/20261025_read_rescue/`;
  2. P1: BAM/VCF headers in v1 contain local absolute paths (e.g. 6 in the S04 locus BAM, 4 in the Manta PASS VCF, 7 in the chr11 SNV VCF);
  3. P2–P6: added client-request coverage, key-event read ledger, primer design, data inventory, external-input status, §9 note.
- **How to compare:** `diff <(cut -f1-8 v1/.../tables/junction_table_target.tsv) <(cut -f1-8 ../../../custom_research_report_supplementary_20261025/tables/junction_table_target.tsv)`;
  report: `diff _client_report_source_v1.md ../../../docs/_client_report_source.md`.
- Scripts that produced v1: git commit `9df28f8` (projects_2026Q4 repo).

- **Note (2026-10-25):** `custom_research_report_supplementary_20261025/` here is the `cp -a` copy; `..._original_moved/` is the
  original folder moved out of the project root so the round-2 build starts from an empty folder (no stale files).
  The two were verified identical with `diff -r` before the move; one of them can be removed with the user's consent.
