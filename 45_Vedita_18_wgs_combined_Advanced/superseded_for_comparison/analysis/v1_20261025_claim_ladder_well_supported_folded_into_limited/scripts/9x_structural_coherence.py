"""Stage 9x — claim ladder (P45 plan section 4.2, pre-registered and approved 2026-10-25 as Q2; not changed after results).

Input: inherited P44 round-2 Stage 7 junction_annotated.tsv (evidence-strength class, unchanged) and Stage 6
junction_by_sample.tsv; W-independent Stage 9 depth-step estimates (MAPQ >= 20) and chr11 BAF segments.

claim_level (one per junction; P44 classification is kept as a separate column, never replaced):
  not_supported                    P44 class Not_supported (no fragment after competitive re-genotyping)
  mapping_observation              max after-rescue fragments in any reported sample <= 1
  repeat_ambiguous_junction        P44 class Ambiguous_repeat, except the partner-end-only case below
  limited_junction_candidate       P44 class Limited (or Well_supported) that is not structurally coherent
  structurally_coherent_limited_candidate
                                   P44 class Limited / Well_supported / Ambiguous_repeat with ambiguity_reason ==
                                   'partner_placement' (the W-side flank is unique), AND in at least one reported
                                   sample with >= 2 after-rescue fragments:
                                   (1) a chr11 breakpoint interval of the junction intersects that sample's depth-step
                                       95% CI [ci95_lo, ci95_hi + 200) (200 bp bins; MAPQ >= 20 estimate), and
                                   (2) that sample's distal chr11q segment (114.7-135.0 Mb) is LOH-like at >= 90% of
                                       baseline-heterozygous sites.
  confirmed                        requires an independent assay (PCR / FISH raw results). None was provided, so this
                                   level is never assigned by this script.
Output: analysis/20261025_p45/stage9x_coherence/claim_ladder.tsv (+ per-sample test rows: coherence_tests.tsv)
"""
import collections
import os
import sys
sys.path.insert(0, os.path.dirname(__file__))
import p45lib as L

OUT = f"{L.OUT}/stage9x_coherence"
os.makedirs(OUT, exist_ok=True)
REP = [r["sample_id"] for r in L.samples(reported_only=True)]
LOH_MIN = 0.90
BIN = 200


def tsv(path):
    with open(path) as fh:
        h = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            yield dict(zip(h, line.rstrip("\n").split("\t")))


ann = list(tsv(f"{L.P44R2}/stage7_annotation/junction_annotated.tsv"))
by = collections.defaultdict(dict)
for d in tsv(f"{L.P44R2}/stage6_genotype/junction_by_sample.tsv"):
    by[d["junction_id"]][d["sample"]] = d
step = {d["sample"]: d for d in tsv(f"{L.SUP}/depth_boundary/boundary_estimates.tsv") if d["mapq_min"] == "20"}
loh = {d["sample"]: float(d["frac_LOH_like"]) for d in tsv(f"{L.SUP}/chr11_baf/baf_segments.tsv")
       if d["segment"] == "distal_114.7-135.0Mb"}


def chr11_intervals(a):
    out = []
    for side in ("A", "B"):
        c = a[f"breakend_{side}"].rsplit(":", 2)[0]
        if c == L.CHR:
            lo, hi = map(int, a[f"uncertainty_{side}"].split("-"))
            out.append((side, lo, hi))
    return out


rows, tests = [], []
for a in ann:
    jid, cls = a["junction_id"], a["classification"]
    maxf = int(a["max_after_frag_S01_S10"])
    eligible = cls in ("Limited", "Well_supported") or (cls == "Ambiguous_repeat" and a["ambiguity_reason"] == "partner_placement")
    coherent_in = []
    if eligible and maxf >= 2:
        for s in REP:
            if s not in by[jid] or int(by[jid][s]["after_frag"]) < 2:
                continue
            st = step.get(s)
            ci = (int(st["ci95_lo"]), int(st["ci95_hi"]) + BIN) if st else None
            hit = [f"{side}:{lo}-{hi}" for side, lo, hi in chr11_intervals(a) if ci and lo < ci[1] and hi >= ci[0]]
            c1 = bool(hit)
            c2 = loh.get(s, 0.0) >= LOH_MIN
            tests.append((jid, s, by[jid][s]["after_frag"], f"{ci[0]}-{ci[1] - 1}" if ci else "NA",
                          ";".join(hit), int(c1), loh.get(s, "NA"), int(c2), int(c1 and c2)))
            if c1 and c2:
                coherent_in.append(s.split("_")[0])
    if cls == "Not_supported":
        lvl = "not_supported"
    elif maxf <= 1:
        lvl = "mapping_observation"
    elif coherent_in:
        lvl = "structurally_coherent_limited_candidate"
    elif cls == "Ambiguous_repeat":
        lvl = "repeat_ambiguous_junction"
    else:
        lvl = "limited_junction_candidate"
    rows.append((jid, cls, a["ambiguity_reason"], maxf, lvl, ",".join(coherent_in)))

with open(f"{OUT}/claim_ladder.tsv", "w") as fh:
    fh.write("junction_id\tevidence_class\tambiguity_reason\tmax_after_frag_S01_S10\tclaim_level\tcoherent_in_samples\n")
    for r in rows:
        fh.write("\t".join(map(str, r)) + "\n")
with open(f"{OUT}/coherence_tests.tsv", "w") as fh:
    fh.write("junction_id\tsample\tafter_frag\tdepth_step_ci95\tchr11_intervals_in_ci\tcriterion1_in_step_ci\t"
             "distal_frac_LOH_like\tcriterion2_LOH_ge_0.90\tcoherent\n")
    for r in tests:
        fh.write("\t".join(map(str, r)) + "\n")
assert not any(r[4] == "confirmed" for r in rows)
print("stage9x OK", collections.Counter(r[4] for r in rows))
for r in rows:
    if r[4] not in ("not_supported", "mapping_observation"):
        print("  ", r)
