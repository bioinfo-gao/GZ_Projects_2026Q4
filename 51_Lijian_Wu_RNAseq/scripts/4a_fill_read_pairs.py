#!/usr/bin/env python3
"""P51 — fill scripts/0_sample_table.csv `raw_read_pairs` from nf-core FastQC (raw, R1).
Run: /Work_bio/gao/configs/.conda/envs/regular_bioinfo/bin/python scripts/4a_fill_read_pairs.py
Source = 'Total Sequences' in output_results/fastqc/raw/<sample>_1_fastqc.zip (R1 of the pair;
R2 must match, asserted). Exits non-zero if any of the 24 samples cannot be filled."""
import csv, glob, os, re, sys, zipfile
S = os.path.dirname(os.path.abspath(__file__)); P = os.path.dirname(S)
T = os.path.join(S, "0_sample_table.csv")

def total(z):
    with zipfile.ZipFile(z) as zf:
        name = next(n for n in zf.namelist() if n.endswith("fastqc_data.txt"))
        for line in zf.read(name).decode().splitlines():
            if line.startswith("Total Sequences"):
                return int(line.split("\t")[1])
    raise ValueError(z)

rows = list(csv.DictReader(open(T)))
missing = []
for d in rows:
    sid = d["sample_id"]
    z1 = glob.glob(os.path.join(P, "output_results/fastqc/raw", f"{sid}_raw_1_fastqc.zip")) or \
         glob.glob(os.path.join(P, "output_results/fastqc/raw", f"{sid}_1_fastqc.zip"))
    z2 = glob.glob(os.path.join(P, "output_results/fastqc/raw", f"{sid}_raw_2_fastqc.zip")) or \
         glob.glob(os.path.join(P, "output_results/fastqc/raw", f"{sid}_2_fastqc.zip"))
    if not (z1 and z2):
        missing.append(sid); continue
    n1, n2 = total(z1[0]), total(z2[0])
    assert n1 == n2, f"{sid}: R1 {n1} != R2 {n2}"
    d["raw_read_pairs"] = str(n1)
with open(T, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
tot = sum(int(d["raw_read_pairs"]) for d in rows if d["raw_read_pairs"])
print(f"filled {len(rows)-len(missing)}/{len(rows)}; total read pairs {tot/1e6:.1f} M; missing {missing}")
sys.exit(1 if missing else 0)
