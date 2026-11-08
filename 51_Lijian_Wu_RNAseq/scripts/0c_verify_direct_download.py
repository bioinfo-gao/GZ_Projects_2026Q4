#!/usr/bin/env python3
"""
P51 — verify FASTQs fetched by 0b_direct_download_dropbox.sh against Dropbox's own content_hash
(SHA-256 over the concatenated SHA-256 digests of each 4 MiB block; Dropbox API spec).
Run: /Work_bio/gao/configs/.conda/envs/regular_bioinfo/bin/python scripts/0c_verify_direct_download.py
Input : logs/direct_download_expected_20261108.tsv  (name, size, content_hash — from the Dropbox API)
Output: logs/direct_download_verified.tsv            (name, size, content_hash, status, checked_at)
Only rows with status=ok are used by 0_build_sample_table.py. Already-ok files are not re-hashed.
"""
import csv, datetime as dt, hashlib, os, sys
from concurrent.futures import ThreadPoolExecutor
PROJ = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DIR = os.path.join(PROJ, "raw_fastq_direct_dropbox_20261108")
EXP = os.path.join(PROJ, "logs", "direct_download_expected_20261108.tsv")
OUT = os.path.join(PROJ, "logs", "direct_download_verified.tsv")
BLOCK = 4 * 1024 * 1024

def dbx_hash(path):
    outer = hashlib.sha256()
    with open(path, "rb") as fh:
        while True:
            b = fh.read(BLOCK)
            if not b: break
            outer.update(hashlib.sha256(b).digest())
    return outer.hexdigest()

done = {}
if os.path.exists(OUT):
    with open(OUT) as fh:
        for r in csv.reader(fh, delimiter="\t"):
            done[r[0]] = r
exp = [r for r in csv.reader(open(EXP), delimiter="\t")]

def check(row):
    name, size, h = row
    p = os.path.join(DIR, name)
    if name in done and done[name][3] == "ok": return done[name]
    if not os.path.exists(p): return [name, size, h, "absent", ""]
    if str(os.path.getsize(p)) != size: return [name, size, h, "size_mismatch", ""]
    st = "ok" if dbx_hash(p) == h else "HASH_MISMATCH"
    return [name, size, h, st, dt.datetime.now().strftime("%F %T")]

with ThreadPoolExecutor(max_workers=4) as ex:
    res = list(ex.map(check, exp))
with open(OUT, "w", newline="") as fh:
    csv.writer(fh, delimiter="\t").writerows(res)
from collections import Counter
c = Counter(r[3] for r in res)
print(f"verified: {dict(c)}  ({len(res)} expected)")
for r in res:
    if r[3] == "HASH_MISMATCH": print("  !! HASH_MISMATCH", r[0])
sys.exit(1 if c.get("HASH_MISMATCH") else 0)
