# v3_20261025 — END taken from rec.stop (still wrong for some BCFs)

- Produced 2026-10-25 16:14–16:20. Original files, moved.
- Defect: v3 used rec.stop when END is absent from rec.info. rec.stop comes from the BCF record's rlen field, which DELLY
  writes as 1 for symbolic alleles in its genotyped BCFs (e.g. DEL00003167 END=157172208 but rec.stop = POS). The P41
  discovery BCFs carry the correct rlen, so v3 was right for them and wrong for re-genotyped DELLY files.
- Fix (v4): END is parsed from the record's INFO text (END=...), falling back to rec.stop only if absent.
