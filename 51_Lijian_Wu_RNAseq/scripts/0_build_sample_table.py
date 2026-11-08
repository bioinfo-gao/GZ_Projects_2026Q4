#!/usr/bin/env python3
"""
P51 — Lijian Wu mouse bulk RNA-seq (control vs Crtc2 line x un/TNF/PGE2/adenosine, n=3).
Build scripts/0_sample_table.csv (SSOT) from the client's registration form + the FASTQ
files actually on disk.

Run:  /Work_bio/gao/configs/.conda/envs/regular_bioinfo/bin/python scripts/0_build_sample_table.py
      (do NOT use `conda run -n regular_bioinfo`: two envs share that name on this machine and
       `-n` resolves to one without openpyxl — software_pitfalls/conda_env_duplicate_names.md)

Re-run any time: Dropbox syncs progressively, so file presence / size is re-measured each run
and stamped into `measured_at`. A sample is `paired_on_disk=yes` only when BOTH mates exist.

Client name -> FASTQ token mapping (assumption, recorded in docs/P51_research_plan_and_review_1107.md §2):
  form Section 4 "control-un-N"        <-> FASTQ "Control_N"
  form "control-PEG2-N"                <-> FASTQ "Control_PGE_N"   (PEG2 = typo of PGE2;
                                            form Section 3G itself writes "pge2")
  form "control-TNF-N"                 <-> FASTQ "Control_TNFA_N"
  form "control-adenosine-N"           <-> FASTQ "Control_Ad_N"
  form "Crtc2-*"                       <-> FASTQ "Crtc2_*" (same pattern)
  replicate N in the form <-> replicate N in the FASTQ name.
"""
import csv
import datetime as dt
import glob
import os
import re
import sys

import openpyxl

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJ_DIR = os.path.dirname(SCRIPT_DIR)
FORM = os.path.join(PROJ_DIR, "docs", "Sample_Registration_Form.xlsx")
FASTQ_DIR = "/Work_bio/dropbox/Dropbox_Data/Quote_09292601"
OUT = os.path.join(SCRIPT_DIR, "0_sample_table.csv")

# treatment token in the client's sample names -> (analysis label, FASTQ token)
TREAT = {
    "un":        ("un",        ""),
    "PEG2":      ("PGE2",      "_PGE"),
    "TNF":       ("TNF",       "_TNFA"),
    "adenosine": ("adenosine", "_Ad"),
}
GENO = {"control": ("control", "Control"), "Crtc2": ("Crtc2", "Crtc2")}

wb = openpyxl.load_workbook(FORM, data_only=True)
ws = wb["Sample Registration"]

# Section 4 per-sample rows: header row has 'Sample Name*' in column A
hdr_row = next(r for r in range(1, ws.max_row + 1) if ws.cell(r, 1).value == "Sample Name*")
form_rows = []
r = hdr_row + 1
while ws.cell(r, 1).value not in (None, ""):
    form_rows.append([ws.cell(r, c).value for c in range(1, 10)])
    r += 1

now = dt.datetime.now().strftime("%Y-%m-%d %H:%M")
rows, problems = [], []
for v in form_rows:
    client_name = str(v[0]).strip()
    m = re.fullmatch(r"(control|Crtc2)-(un|PEG2|TNF|adenosine)-([123])", client_name)
    if not m:
        problems.append(f"unparsed client sample name: {client_name!r}")
        continue
    geno_c, treat_c, rep = m.groups()
    geno, geno_tok = GENO[geno_c]
    treat, treat_tok = TREAT[treat_c]
    group = f"{geno}-{treat}"
    sample_id = f"{group}-{rep}"
    token = f"{geno_tok}{treat_tok}_{rep}"
    r1 = glob.glob(os.path.join(FASTQ_DIR, f"*_{token}_S*_L*_R1_001.fastq.gz"))
    r2 = glob.glob(os.path.join(FASTQ_DIR, f"*_{token}_S*_L*_R2_001.fastq.gz"))
    if len(r1) > 1 or len(r2) > 1:
        problems.append(f"{sample_id}: >1 file per mate for token {token} — multi-lane not expected")
    r1 = r1[0] if r1 else ""
    r2 = r2[0] if r2 else ""
    snum = re.search(r"_S(\d+)_L", r1 or r2)
    gib = lambda f: f"{os.path.getsize(f) / 2**30:.2f}" if f else ""
    rows.append({
        "sample_id": sample_id,
        "client_sample_name": client_name,
        "group": group,
        "genotype": geno,
        "treatment": treat,
        "replicate": rep,
        "fastq_token": token,
        "illumina_S": f"S{snum.group(1)}" if snum else "",
        "species": "Mus musculus (GRCm39)",
        "sample_material": "Cell line (established); Normal only; Snap-frozen (per form 3C)",
        "service": v[1],
        "paired_on_disk": "yes" if (r1 and r2) else "no",
        "r1_file": os.path.basename(r1),
        "r2_file": os.path.basename(r2),
        "R1_GiB": gib(r1),
        "R2_GiB": gib(r2),
        "raw_read_pairs": "",          # filled from pipeline FastQC after the run
        "fastq_dir": FASTQ_DIR,
        "measured_at": now,
    })

# keep the read-pair column if a previous run already filled it
if os.path.exists(OUT):
    with open(OUT) as fh:
        old = {d["sample_id"]: d for d in csv.DictReader(fh)}
    for d in rows:
        if d["sample_id"] in old and old[d["sample_id"]].get("raw_read_pairs"):
            if old[d["sample_id"]].get("r1_file") == d["r1_file"]:
                d["raw_read_pairs"] = old[d["sample_id"]]["raw_read_pairs"]

with open(OUT, "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
    w.writeheader()
    w.writerows(rows)

n_pair = sum(d["paired_on_disk"] == "yes" for d in rows)
tot = sum(float(d["R1_GiB"] or 0) + float(d["R2_GiB"] or 0) for d in rows)
print(f"[{now}] form samples: {len(rows)} | paired on disk: {n_pair}/{len(rows)} | "
      f"FASTQ on disk: {tot:.1f} GiB")
for d in rows:
    flag = "" if d["paired_on_disk"] == "yes" else \
        f"  <-- missing {'R1 ' if not d['r1_file'] else ''}{'R2' if not d['r2_file'] else ''}"
    print(f"  {d['sample_id']:<22} {d['client_sample_name']:<22} {d['illumina_S']:<4} "
          f"{d['R1_GiB']:>6} {d['R2_GiB']:>6}{flag}")

# every FASTQ in the folder must be claimed by exactly one form sample
claimed = {d["r1_file"] for d in rows} | {d["r2_file"] for d in rows}
orphans = [os.path.basename(f) for f in glob.glob(os.path.join(FASTQ_DIR, "*.fastq.gz"))
           if os.path.basename(f) not in claimed]
if orphans:
    problems.append(f"{len(orphans)} FASTQ not matched to any form sample: {orphans}")
assert len(rows) == 24, f"expected 24 form samples, got {len(rows)}"
assert len({d['sample_id'] for d in rows}) == 24
assert all(sum(x["group"] == g for x in rows) == 3 for g in {d["group"] for d in rows}), "not 8x3"
if problems:
    print("\n!! PROBLEMS:"); [print("  -", p) for p in problems]
    sys.exit(1)
print(f"OK -> {OUT}")
