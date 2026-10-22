#!/usr/bin/env python3
"""Classify independent re-alignments against every frozen P41 candidate breakpoint."""
import csv, pathlib, sys
import pysam

candidates_tsv, manifest_tsv, sam_path, outdir = map(pathlib.Path, sys.argv[1:5])
outdir.mkdir(parents=True, exist_ok=True)
candidates = {r["candidate_id"]: r for r in csv.DictReader(candidates_tsv.open(), delimiter="\t")}
manifest = list(csv.DictReader(manifest_tsv.open(), delimiter="\t"))
records = {r["fastq_id"]: [] for r in manifest}
with pysam.AlignmentFile(sam_path) as sam:
    for read in sam:
        if read.query_name in records and not read.is_unmapped:
            records[read.query_name].append(read)

def end_match(read, chrom, pos, direction):
    if read.reference_name != chrom: return False
    end = read.reference_end if direction == "+" else read.reference_start + 1
    return abs(end - int(pos)) <= 20

by_candidate = {}
for item in manifest:
    candidate = candidates[item["candidate_id"]]
    aligns = records[item["fastq_id"]]
    left = [a for a in aligns if end_match(a, candidate["chr1"], candidate["pos1"], candidate["dir1"])]
    right = [a for a in aligns if end_match(a, candidate["chr2"], candidate["pos2"], candidate["dir2"])]
    pairs = [(a, b) for a in left for b in right if a is not b]
    if pairs:
        best = max(pairs, key=lambda z: (z[0].get_tag("AS") if z[0].has_tag("AS") else 0) + (z[1].get_tag("AS") if z[1].has_tag("AS") else 0))
        unique = best[0].mapping_quality >= 20 and best[1].mapping_quality >= 20
        call = "unique_reconstruction" if unique else "ambiguous_reconstruction"
    else:
        call = "not_reconstructed"
    by_candidate.setdefault(item["candidate_id"], []).append((item, call, len(aligns)))

rows = []
for cid, c in candidates.items():
    calls = by_candidate.get(cid, [])
    fragments = {x[0]["qname"] for x in calls}
    unique = len({x[0]["qname"] for x in calls if x[1] == "unique_reconstruction"})
    ambiguous = len({x[0]["qname"] for x in calls if x[1] == "ambiguous_reconstruction"})
    observed = len(fragments)
    if observed == 0: status = "not_evaluable"
    elif unique >= 2: status = "multi_fragment_supported"
    elif unique == 1: status = "single_fragment_supported"
    elif ambiguous: status = "ambiguous_reconstruction"
    else: status = "not_reconstructed"
    rows.append({**c, "source_fragments": c["n_split_all"], "reextracted_fragments": observed,
                 "unique_reconstructed_fragments": unique, "ambiguous_reconstructed_fragments": ambiguous,
                 "rescue_status": status})
fields = list(rows[0])
with (outdir / "rescue_candidate_results.tsv").open("w", newline="") as h:
    w = csv.DictWriter(h, fieldnames=fields, delimiter="\t", lineterminator="\n"); w.writeheader(); w.writerows(rows)
with (outdir / "rescue_sample_summary.tsv").open("w", newline="") as h:
    w = csv.writer(h, delimiter="\t", lineterminator="\n"); w.writerow(["sample", "candidates", "multi_fragment_supported", "single_fragment_supported", "ambiguous_reconstruction", "not_reconstructed", "not_evaluable"])
    for sample in sorted({r["sample"] for r in rows}):
        rr = [r for r in rows if r["sample"] == sample]
        w.writerow([sample, len(rr), *[sum(r["rescue_status"] == x for r in rr) for x in ("multi_fragment_supported", "single_fragment_supported", "ambiguous_reconstruction", "not_reconstructed", "not_evaluable")]])
print(f"classified {len(rows)} candidates")
