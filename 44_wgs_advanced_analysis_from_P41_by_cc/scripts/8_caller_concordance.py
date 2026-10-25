"""Stage 8 — Manta / DELLY concordance for every junction (plan section 3), all 18 samples, from the P41 call sets.

Call sets (read-only, P41):
  manta_pass      per-sample genome-wide diploidSV.vcf.gz, FILTER=PASS          (filtered)
  manta_all       per-sample diploidSV.vcf.gz (any FILTER) + candidateSV.vcf.gz  (unfiltered)
  manta_joint     joint 18-sample target-window diploidSV (graphNodeMaxEdgeCount=1000), sample GT non-ref
  delly_pass      per-sample delly discovery BCF, FILTER=PASS                    (filtered)
  delly_all       per-sample delly discovery BCF (any FILTER)                     (unfiltered)
  delly_regeno    window re-genotyping merged BCF, sample GT non-ref
Each record is converted to a breakend pair (VCF BND bracket notation; Manta SVTYPE; DELLY CT). Match classes:
  two_end    both breakends within TOL bp of the junction's two breakends and orientations equal
  two_end_pos both within TOL but orientation differs
  one_end    one breakend within TOL (anchor-only; not concordance)
  none
Point insertions (SVTYPE=INS) match only junctions with both breakends within TOL of the insertion position.
Output: analysis/<RUN>/stage8_concordance{L.SUFFIX}/{concordance_long.tsv, concordance_summary.tsv}
"""
import collections
import gzip
import os
import re
import sys
import pysam
sys.path.insert(0, os.path.dirname(__file__))
import p44lib as L

RUN = os.environ.get("P44_RUN_DATE", "20261024")
BASE = f"{L.PROJ}/analysis/{RUN}_read_rescue"
OUT = f"{BASE}/stage8_concordance{L.SUFFIX}"
os.makedirs(OUT, exist_ok=True)
P41 = "/home/gao/projects_2026M7/41_Vedict_18_wgs/analysis/sv"
TOL = 100
SAMPLES = [r["sample_id"] for r in L.samples()]
pysam.set_verbosity(0)


def tsv(path):
    with open(path) as fh:
        h = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            yield dict(zip(h, line.rstrip("\n").split("\t")))


def bk(s):
    c, p, o = s.rsplit(":", 2)
    return (c, int(p), o)


BND_RE = re.compile(r"([\[\]])([^:\[\]]+):(\d+)([\[\]])")


def pairs(rec):
    """breakend pairs [(A,B)] or ('INS', chrom, pos)"""
    hdr = rec.header.info

    class _I:                                   # header-safe INFO access (pysam raises on undeclared keys)
        def get(self, k, d=None):
            return rec.info[k] if (k in hdr and k in rec.info) else d

        def __contains__(self, k):
            return k in hdr and k in rec.info
    info = _I()
    t = info.get("SVTYPE", "")
    c, p = rec.chrom, rec.pos
    end = info.get("END", p)
    ct = info.get("CT")
    if t == "INS":
        return [("INS", c, p)]
    if t == "BND":
        alt = str(rec.alts[0]) if rec.alts else ""
        m = BND_RE.search(alt)
        if m:
            c2, p2 = m.group(2), int(m.group(3))
            if alt.startswith(("[", "]")):          # ]p]t or [p[t : piece joined BEFORE t
                a = (c, p, "-")
                b = (c2, p2, "+" if m.group(1) == "]" else "-")
            else:                                     # t[p[ or t]p]
                a = (c, p, "+")
                b = (c2, p2, "-" if m.group(1) == "[" else "+")
            return [(a, b)]
        if ct and info.get("CHR2"):                  # DELLY BND
            c2, p2 = info["CHR2"], info.get("POS2", info.get("END"))
            oa, ob = {"3to5": ("+", "-"), "5to3": ("-", "+"), "3to3": ("+", "+"), "5to5": ("-", "-")}[ct]
            return [((c, p, oa), (c2, p2, ob))]
        return []
    if ct:                                            # DELLY DEL/DUP/INV
        oa, ob = {"3to5": ("+", "-"), "5to3": ("-", "+"), "3to3": ("+", "+"), "5to5": ("-", "-")}[ct]
        return [((c, p, oa), (c, end, ob))]
    if t == "DEL":
        return [((c, p, "+"), (c, end + 1, "-"))]
    if t == "DUP":
        return [((c, p, "-"), (c, end, "+"))]
    if t == "INV":
        if "INV3" in info:
            return [((c, p, "+"), (c, end, "+"))]
        if "INV5" in info:
            return [((c, p + 1, "-"), (c, end + 1, "-"))]
        return [((c, p, "+"), (c, end, "+")), ((c, p + 1, "-"), (c, end + 1, "-"))]
    return []


def classify(J, recpairs):
    A, B = J
    best = "none"
    rank = {"none": 0, "one_end": 1, "two_end_pos": 2, "two_end": 3}
    for pr in recpairs:
        if pr[0] == "INS":
            _, c, p = pr
            cls = "two_end" if (A[0] == B[0] == c and abs(A[1] - p) <= TOL and abs(B[1] - p) <= TOL) else \
                  ("one_end" if any(x[0] == c and abs(x[1] - p) <= TOL for x in (A, B)) else "none")
        else:
            a, b = pr
            cls = "none"
            for x, y in ((a, b), (b, a)):
                if x[0] == A[0] and y[0] == B[0] and abs(x[1] - A[1]) <= TOL and abs(y[1] - B[1]) <= TOL:
                    cls = "two_end" if (x[2] == A[2] and y[2] == B[2]) else "two_end_pos"
                    break
            if cls == "none" and any(z[0] == w[0] and abs(z[1] - w[1]) <= TOL for z in (a, b) for w in (A, B)):
                cls = "one_end"
        if rank[cls] > rank[best]:
            best = cls
    return best


def fetch(vf, b, sample=None, pass_only=False):
    out = []
    try:
        it = vf.fetch(b[0], max(0, b[1] - TOL - 1), b[1] + TOL)
    except (ValueError, KeyError):
        return out
    for rec in it:
        if pass_only and "PASS" not in rec.filter.keys():
            continue
        if sample is not None:
            gt = rec.samples[sample].get("GT")
            if gt is None or all(g in (0, None) for g in gt):
                continue
        out.append((rec, pairs(rec)))
    return out


callsets = {}
for s in SAMPLES:
    v = f"{P41}/manta_gw/{s}/results/variants"
    callsets[(s, "manta_pass")] = [(f"{v}/diploidSV.vcf.gz", None, True)]
    callsets[(s, "manta_all")] = [(f"{v}/diploidSV.vcf.gz", None, False), (f"{v}/candidateSV.vcf.gz", None, False)]
    callsets[(s, "manta_joint")] = [(f"{P41}/manta_region_joint_20261007/results/variants/diploidSV.vcf.gz", s, False)]
    callsets[(s, "delly_pass")] = [(f"{P41}/delly/discover/{s}.delly.bcf", None, True)]
    callsets[(s, "delly_all")] = [(f"{P41}/delly/discover/{s}.delly.bcf", None, False)]
    callsets[(s, "delly_regeno")] = [(f"{P41}/delly/region_geno_20261007/merged_genotypes.bcf", s, False)]
handles = {}


def vfile(path):
    if path not in handles:
        handles[path] = pysam.VariantFile(path)
    return handles[path]


junctions = list(tsv(f"{BASE}/stage7_annotation{L.SUFFIX}/junction_annotated.tsv"))
long_rows = []
summary = []
for j in junctions:
    A, B = bk(j["breakend_A"]), bk(j["breakend_B"])
    per = collections.defaultdict(dict)
    for (s, cs), files in callsets.items():
        recs = []
        ids = []
        for path, smp, pss in files:
            for b in (A, B):
                for rec, pr in fetch(vfile(path), b, smp, pss):
                    recs += pr
                    ids.append(rec.id or f"{rec.chrom}:{rec.pos}")
        cls = classify((A, B), recs)
        per[s][cs] = cls
        long_rows.append((j["junction_id"], s, cs, cls, ";".join(sorted(set(ids))[:5])))
    rep = [s for s in SAMPLES if int(s[1:3]) <= 10]
    m = lambda cs, ss: [x.split("_")[0] for x in ss if per[x][cs] == "two_end"]
    manta_any = m("manta_all", rep) + m("manta_joint", rep)
    delly_any = m("delly_all", rep) + m("delly_regeno", rep)
    callers = "Manta+DELLY" if manta_any and delly_any else "Manta" if manta_any else "DELLY" if delly_any else "neither"
    summary.append(dict(junction_id=j["junction_id"], caller_support=callers,
                        manta_pass_samples=",".join(m("manta_pass", rep)),
                        manta_unfiltered_or_joint_samples=",".join(sorted(set(manta_any))),
                        delly_pass_samples=",".join(m("delly_pass", rep)),
                        delly_unfiltered_or_regeno_samples=",".join(sorted(set(delly_any))),
                        any_one_end_anchor=",".join(sorted({x.split("_")[0] for x in rep for cs in per[x]
                                                            if per[x][cs] in ("one_end", "two_end_pos")}))))
with open(f"{OUT}/concordance_long.tsv", "w") as fh:
    fh.write("junction_id\tsample\tcallset\tmatch\trecord_ids\n")
    for r in long_rows:
        fh.write("\t".join(r) + "\n")
with open(f"{OUT}/concordance_summary.tsv", "w") as fh:
    fh.write("\t".join(summary[0].keys()) + "\n")
    for r in summary:
        fh.write("\t".join(r.values()) + "\n")
print("stage8 OK", collections.Counter(r["caller_support"] for r in summary))
