#!/usr/bin/env python3
"""Build the P42 candidate SSOT from frozen P41 target-evidence tables."""
import csv, pathlib, sys

samples_csv, evidence_dir, outdir = map(pathlib.Path, sys.argv[1:4])
outdir.mkdir(parents=True, exist_ok=True)
samples = list(csv.DictReader(samples_csv.open()))
fields = ["candidate_id", "sample", "type", "chr1", "pos1", "dir1", "chr2", "pos2", "dir2",
          "n_split_all", "n_split_nondup", "n_split_nondup_mapq20", "source_reads"]
rows, support = [], []
for sample in samples:
    sid = sample["sample_id"]
    source = evidence_dir / f"{sid}.junctions.tsv"
    if not source.is_file(): raise SystemExit(f"missing P41 source: {source}")
    eligible = []
    with source.open() as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["in_target"] == "True" and int(row["n_split_nondup_mapq20"]) >= 1:
                eligible.append(row)
    for i, row in enumerate(eligible, 1):
        cid = f"{sid}_J{i:03d}"
        qnames = [x for x in row["reads"].split(",") if x]
        if len(qnames) != int(row["n_split_all"]):
            raise SystemExit(f"{cid}: source qname list is incomplete ({len(qnames)} != {row['n_split_all']})")
        rows.append(dict(candidate_id=cid, sample=sid, type=row["type"], chr1=row["chr1"], pos1=row["pos1"],
                         dir1=row["dir1"], chr2=row["chr2"], pos2=row["pos2"], dir2=row["dir2"],
                         n_split_all=row["n_split_all"], n_split_nondup=row["n_split_nondup"],
                         n_split_nondup_mapq20=row["n_split_nondup_mapq20"], source_reads=row["reads"]))
        support.extend({"candidate_id": cid, "sample": sid, "qname": q} for q in qnames)
    expected = int(sample["source_candidate_count"])
    if len(eligible) != expected:
        raise SystemExit(f"{sid}: expected {expected} eligible candidates, found {len(eligible)}")
for required in ("S02_DLD1_c1_ctrl", "S03_DLD1_c2_ctrl", "S10_GM12878_c8"):
    if not any(r["sample"] == required for r in rows): raise SystemExit(f"required sample absent: {required}")
if len(rows) != 179: raise SystemExit(f"candidate total {len(rows)} != 179")
with (outdir / "candidate_junctions.tsv").open("w", newline="") as h:
    w = csv.DictWriter(h, fieldnames=fields, delimiter="\t", lineterminator="\n"); w.writeheader(); w.writerows(rows)
with (outdir / "candidate_support_fragments.tsv").open("w", newline="") as h:
    w = csv.DictWriter(h, fieldnames=["candidate_id", "sample", "qname"], delimiter="\t", lineterminator="\n"); w.writeheader(); w.writerows(support)
with (outdir / "candidate_counts.tsv").open("w", newline="") as h:
    w = csv.writer(h, delimiter="\t", lineterminator="\n"); w.writerow(["sample", "candidate_count"])
    for s in samples: w.writerow([s["sample_id"], sum(r["sample"] == s["sample_id"] for r in rows)])
print(f"prepared {len(rows)} candidates and {len(support)} source-fragment links")
