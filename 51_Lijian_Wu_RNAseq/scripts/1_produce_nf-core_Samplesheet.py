#!/usr/bin/env python3
"""
P51 — build nf-core/rnaseq samplesheet from scripts/0_sample_table.csv (SSOT).

Run:  /Work_bio/gao/configs/.conda/envs/regular_bioinfo/bin/python scripts/1_produce_nf-core_Samplesheet.py [--require-all]

- Only samples with BOTH mates on disk AND passing `gzip -t` are emitted (Dropbox syncs
  progressively; a partially-present sample must never reach nextflow).
- gzip results are cached in logs/gzip_check.tsv keyed by (file, size, mtime): re-runs only
  test files that are new or changed.
- --require-all: exit 1 unless all 24 designed samples are emitted (used by the auto-launcher
  before the full run).
- strandedness = auto (nf-core infers via Salmon subsample; client form does not state it).
- Output: scripts/nf_core_samplesheet.csv ; last line "EMITTED <n> <total>" for the launcher.
"""
import csv
import os
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJ_DIR = os.path.dirname(SCRIPT_DIR)
TABLE = os.path.join(SCRIPT_DIR, "0_sample_table.csv")
OUT = os.path.join(SCRIPT_DIR, "nf_core_samplesheet.csv")
CACHE = os.path.join(PROJ_DIR, "logs", "gzip_check.tsv")
REQUIRE_ALL = "--require-all" in sys.argv

os.makedirs(os.path.dirname(CACHE), exist_ok=True)
cache = {}
if os.path.exists(CACHE):
    with open(CACHE) as fh:
        for line in fh:
            f, size, mtime, ok = line.rstrip("\n").split("\t")
            cache[(f, size, mtime)] = ok

with open(TABLE) as fh:
    samples = list(csv.DictReader(fh))


def key(path):
    st = os.stat(path)
    return (path, str(st.st_size), str(int(st.st_mtime)))


def gz_ok(path):
    k = key(path)
    if k in cache:
        return k, cache[k]
    rc = subprocess.run(["gzip", "-t", path], stdout=subprocess.DEVNULL,
                        stderr=subprocess.DEVNULL).returncode
    return k, "ok" if rc == 0 else "FAIL"


paths = []
for s in samples:
    if s["paired_on_disk"] == "yes":
        for col in ("r1_file", "r2_file"):
            paths.append(s[col.replace("file", "path")])

with ThreadPoolExecutor(max_workers=8) as ex:
    results = dict(ex.map(gz_ok, paths))
cache.update(results)
with open(CACHE, "w") as fh:
    for (f, size, mtime), ok in sorted(cache.items()):
        fh.write(f"{f}\t{size}\t{mtime}\t{ok}\n")
status = {k[0]: v for k, v in results.items()}

rows, skipped = [], []
for s in samples:
    sid = s["sample_id"]
    if s["paired_on_disk"] != "yes":
        skipped.append(f"{sid}: not paired on disk yet")
        continue
    r1 = s["r1_path"]
    r2 = s["r2_path"]
    bad = [os.path.basename(p) for p in (r1, r2) if status.get(p) != "ok"]
    if bad:
        skipped.append(f"{sid}: gzip -t FAILED {bad}")
        continue
    # guard against a silent sample swap: file must carry this sample's own token
    tok = f"_{s['fastq_token']}_S"
    assert tok in os.path.basename(r1) and tok in os.path.basename(r2), sid
    rows.append({"sample": sid, "fastq_1": r1, "fastq_2": r2, "strandedness": "auto"})

with open(OUT, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=["sample", "fastq_1", "fastq_2", "strandedness"])
    w.writeheader()
    w.writerows(rows)

assert len({r["fastq_1"] for r in rows}) == len(rows), "duplicate FASTQ emitted"
print(f"Written {OUT}: {len(rows)}/{len(samples)} samples")
for x in skipped:
    print("  skipped -", x)
print(f"EMITTED {len(rows)} {len(samples)}")
if REQUIRE_ALL and len(rows) != len(samples):
    sys.exit(1)
