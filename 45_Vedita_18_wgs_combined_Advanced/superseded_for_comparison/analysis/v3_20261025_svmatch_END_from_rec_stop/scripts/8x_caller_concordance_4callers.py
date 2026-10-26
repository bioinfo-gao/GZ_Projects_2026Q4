"""Stage 8x — four-caller concordance for every junction (P45 plan section 3/5): Manta, DELLY (P41 call sets, exactly as
P44 Stage 8) plus GRIDSS (P42 joint 18-sample run) and SvABA (P42 DLD1 S01-S05 multisample run, v5), both run by P42 on
chr11:114,554,000-114,675,500 (= W).

Record -> breakend pair conversion and match classes are copied unchanged from P44 8_caller_concordance.py
(from_P44/), so Manta/DELLY classes are recomputed and ASSERTED equal to the inherited P44 Stage 8 output.
Additions:
  * GRIDSS single breakends (ALT like 'G.' / '.G', no partner) -> anchor only ('one_end' at most).
  * best_two_end_dist: over records whose two breakends both lie within TOL, the smallest max(|dA|, |dB|);
    within_50bp flag = best_two_end_dist <= 50 (the tolerance P42 used), so both rules come from one table.
  * per-sample presence: Manta/DELLY as P44; GRIDSS FORMAT VF (BVF for single breakends) > 0; SvABA FORMAT AD > 0.
Call sets:
  gridss_all   joint VCF, any FILTER, sample VF/BVF > 0          (unfiltered)
  gridss_pass  joint VCF, FILTER=PASS, sample VF/BVF > 0         (filtered)
  svaba_all    unfiltered germline + somatic SV VCF, sample AD > 0 (S01-S05 only)
  svaba_pass   filtered germline + somatic SV VCF, sample AD > 0   (S01-S05 only)
Output: analysis/20261025_p45/stage8x_concordance/{concordance_long_4callers.tsv, concordance_summary_4callers.tsv}
"""
import collections
import os
import re
import sys
import pysam
sys.path.insert(0, os.path.dirname(__file__))
import p45lib as L

OUT = f"{L.OUT}/stage8x_concordance"
os.makedirs(OUT, exist_ok=True)
P41 = "/home/gao/projects_2026M7/41_Vedict_18_wgs/analysis/sv"
EXT = f"{L.P42IN}/external_callers"
TOL = 100
SAMPLES = [r["sample_id"] for r in L.samples()]
REP = [s for s in SAMPLES if int(s[1:3]) <= 10]
DLD1 = [s for s in SAMPLES if int(s[1:3]) <= 5]
pysam.set_verbosity(0)


def tsv(path):
    with open(path) as fh:
        h = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            yield dict(zip(h, line.rstrip("\n").split("\t")))


def bk(s):
    c, p, o = s.rsplit(":", 2)
    return (c, int(p), o)


import svmatch as SM
from svmatch import pairs, classify
SM.TOL = TOL


# ---------------------------------------------------------------- record selection (P44 fetch semantics for every call set)
def carrier(rec, kind, sample):
    """P44 rule for multi-sample GT files; P45 rules for GRIDSS (VF, or BVF for single breakends) and SvABA (AD)."""
    if kind == "single":
        return True
    if kind == "svaba":                               # SvABA sample columns are BAM paths
        col = next(c for c in rec.samples.keys() if os.path.basename(c).startswith(sample))
    else:
        col = sample
    f = rec.samples[col]
    if kind == "gt":
        gt = f.get("GT")
        return not (gt is None or all(g in (0, None) for g in gt))
    if kind == "gridss":
        v = f.get("VF")
        if v is None:
            v = f.get("BVF")
        return bool(v)
    if kind == "svaba":
        v = f.get("AD")
        if isinstance(v, tuple):
            v = v[-1]
        return bool(v)
    raise ValueError(kind)


callsets = {}
for s in SAMPLES:
    v = f"{P41}/manta_gw/{s}/results/variants"
    callsets[(s, "manta_pass")] = [(f"{v}/diploidSV.vcf.gz", "single", True)]
    callsets[(s, "manta_all")] = [(f"{v}/diploidSV.vcf.gz", "single", False), (f"{v}/candidateSV.vcf.gz", "single", False)]
    callsets[(s, "manta_joint")] = [(f"{P41}/manta_region_joint_20261007/results/variants/diploidSV.vcf.gz", "gt", False)]
    callsets[(s, "delly_pass")] = [(f"{P41}/delly/discover/{s}.delly.bcf", "single", True)]
    callsets[(s, "delly_all")] = [(f"{P41}/delly/discover/{s}.delly.bcf", "single", False)]
    callsets[(s, "delly_regeno")] = [(f"{P41}/delly/region_geno_20261007/merged_genotypes.bcf", "gt", False)]
    callsets[(s, "gridss_all")] = [(f"{EXT}/gridss_joint_18samples.vcf.gz", "gridss", False)]
    callsets[(s, "gridss_pass")] = [(f"{EXT}/gridss_joint_18samples.vcf.gz", "gridss", True)]
    if s in DLD1:
        sv = f"{EXT}/svaba_dld1_targeted_multisample.svaba"
        callsets[(s, "svaba_all")] = [(f"{sv}.unfiltered.germline.sv.vcf.gz", "svaba", False),
                                      (f"{sv}.unfiltered.somatic.sv.vcf.gz", "svaba", False)]
        callsets[(s, "svaba_pass")] = [(f"{sv}.germline.sv.vcf.gz", "svaba", False),
                                       (f"{sv}.somatic.sv.vcf.gz", "svaba", False)]

handles = {}


def fetch(path, b, kind, pass_only, sample):
    """records overlapping [pos-TOL, pos+TOL] of breakend b (exactly P44 Stage 8 fetch), filtered by PASS and carrier"""
    if path not in handles:
        handles[path] = pysam.VariantFile(path)
    out = []
    try:
        it = handles[path].fetch(b[0], max(0, b[1] - TOL - 1), b[1] + TOL)
    except (ValueError, KeyError):
        return out
    for rec in it:
        if pass_only and "PASS" not in rec.filter.keys():
            continue
        if not carrier(rec, kind, sample):
            continue
        out.append((rec.id or f"{rec.chrom}:{rec.pos}", pairs(rec)))
    return out


junctions = list(tsv(f"{L.P44R2}/stage7_annotation/junction_annotated.tsv"))
long_rows, summary = [], []
for j in junctions:
    A, B = bk(j["breakend_A"]), bk(j["breakend_B"])
    per = collections.defaultdict(dict)
    for (s, cs), files in callsets.items():
        recs, ids = [], []
        for path, kind, pss in files:
            for b in (A, B):
                for rid, pr in fetch(path, b, kind, pss, s):
                    recs += pr
                    ids.append(rid)
        cls, dist = classify((A, B), recs)
        per[s][cs] = cls
        long_rows.append((j["junction_id"], s, cs, cls, "" if dist is None else str(dist),
                          "1" if (dist is not None and dist <= 50) else "0", ";".join(sorted(set(ids))[:5])))

    def m(cs, ss):
        return [x.split("_")[0] for x in ss if per[x].get(cs) == "two_end"]
    sets = {"Manta": m("manta_all", REP) + m("manta_joint", REP), "DELLY": m("delly_all", REP) + m("delly_regeno", REP),
            "GRIDSS": m("gridss_all", REP), "SvABA": m("svaba_all", DLD1)}
    callers = "+".join(k for k in ("Manta", "DELLY", "GRIDSS", "SvABA") if sets[k]) or "none"
    summary.append(collections.OrderedDict(
        junction_id=j["junction_id"], caller_support_4=callers, n_callers_two_end=sum(1 for k in sets if sets[k]),
        manta_samples=",".join(sorted(set(sets["Manta"]))), manta_PASS_samples=",".join(m("manta_pass", REP)),
        delly_samples=",".join(sorted(set(sets["DELLY"]))), delly_PASS_samples=",".join(m("delly_pass", REP)),
        gridss_samples=",".join(sets["GRIDSS"]), gridss_PASS_samples=",".join(m("gridss_pass", REP)),
        svaba_samples_S01_S05=",".join(sets["SvABA"]), svaba_PASS_samples_S01_S05=",".join(m("svaba_pass", DLD1)),
        any_one_end_anchor=",".join(sorted({x.split("_")[0] for x in REP for cs in per[x]
                                            if per[x][cs] in ("one_end", "two_end_pos")}))))

with open(f"{OUT}/concordance_long_4callers.tsv", "w") as fh:
    fh.write("junction_id\tsample\tcallset\tmatch\tbest_two_end_dist_bp\twithin_50bp\trecord_ids\n")
    for r in long_rows:
        fh.write("\t".join(r) + "\n")
with open(f"{OUT}/concordance_summary_4callers.tsv", "w") as fh:
    fh.write("\t".join(summary[0].keys()) + "\n")
    for r in summary:
        fh.write("\t".join(map(str, r.values())) + "\n")

# ---- reproducibility assertion: Manta/DELLY classes must equal the inherited P44 Stage 8 output exactly
p44 = {(d["junction_id"], d["sample"], d["callset"]): d["match"] for d in tsv(f"{L.P44R2}/stage8_concordance/concordance_long.tsv")}
mine = {(r[0], r[1], r[2]): r[3] for r in long_rows if r[2].split("_")[0] in ("manta", "delly")}
# P45: after the END fix (svmatch.py) differences from P44 Stage 8 are expected for DEL/DUP/INV records; every
# difference is written out, and it must be an upgrade (P44 had END collapsed to POS, so it could only under-match).
assert len(p44) == len(mine), "row sets differ"
diff = [(k, p44[k], mine[k]) for k in p44 if mine.get(k) != p44[k]]
RK = {"none": 0, "one_end": 1, "two_end_pos": 2, "two_end": 3}
with open(f"{OUT}/manta_delly_diff_vs_P44_stage8.tsv", "w") as fh:
    fh.write("junction_id\tsample\tcallset\tP44_stage8_match\tP45_match\n")
    for (j, s, cs), a, b in diff:
        fh.write(f"{j}\t{s}\t{cs}\t{a}\t{b}\n")
down = [d for d in diff if RK[d[2]] < RK[d[1]]]
print(f"Manta/DELLY rows: {len(p44)}; differing from P44 Stage 8 after END fix: {len(diff)}; downgrades: {len(down)}")
print("stage8x OK", collections.Counter(r["caller_support_4"] for r in summary))
