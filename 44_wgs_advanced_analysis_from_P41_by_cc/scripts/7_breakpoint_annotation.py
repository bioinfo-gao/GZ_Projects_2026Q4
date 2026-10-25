"""Stage 7 — base-pair annotation and classification of every precise junction (plan section 4).

Per junction:
  * breakpoint uncertainty from the reference itself: how far the junction can slide left/right while the alt
    sequence stays identical (microhomology / repeat shift); reported as interval at each breakend
  * microhomology sequence, inserted sequence
  * flank uniqueness: 100 bp immediately on each side of the junction, bwa-mem2 -a vs GRCh38; n_equiv = number of
    genomic placements with AS >= best AS - MARGIN (1 = uniquely placeable)
  * RepeatMasker / segmental duplication / simple repeat overlap at each breakend (UCSC hg38 tables)
  * support summary from Stage 6 (S01-S10 reported; S11-S18 background panel)
  * comparison annotation: matching P41 event (both ends within 15 bp, same orientation)
Classification (pre-registered thresholds; plan section 4):
  Well supported   max after-rescue fragments in one sample >= 3; >= 2 evidence types (template split read-ends,
                   discordant pairs, assembled contig); both flanks n_equiv == 1; uncertainty <= 10 bp
  Ambiguous_repeat max fragments >= 2 and (a flank n_equiv > 1 or uncertainty > 10 bp)
  Limited          max fragments >= 2, flanks unique, uncertainty <= 10 bp, but < 3 fragments or 1 evidence type
                   (gap between the plan's Well and Ambiguous levels; added 2026-10-25, see plan change log)
  Weak             max fragments == 1
  Not_supported    0 fragments after competitive re-genotyping in every sample
Output: analysis/<RUN>/stage7_annotation{L.SUFFIX}/{junction_annotated.tsv, flank_alignments.tsv, p41_comparison.tsv}
"""
import collections
import gzip
import os
import subprocess
import sys
import pysam
sys.path.insert(0, os.path.dirname(__file__))
import p44lib as L

RUN = os.environ.get("P44_RUN_DATE", "20261024")
BASE = f"{L.PROJ}/analysis/{RUN}_read_rescue"
OUT = f"{BASE}/stage7_annotation{L.SUFFIX}"
os.makedirs(OUT, exist_ok=True)
ENV = "/Work_bio/gao/configs/.conda/envs/regular_bioinfo/bin"
BWA_IDX = "/Work_bio/references/Homo_sapiens/GRCh38/GATK.GRCh38/BWAmem2Index/Homo_sapiens_assembly38.fasta"
REFS = f"{L.PROJ}/refs"
P41_EVENTS = "/home/gao/projects_2026M7/41_Vedict_18_wgs/analysis/region/compare/events_all.tsv"
FL, MARGIN, MAXSLIDE = 100, 10, 200
fa = pysam.FastaFile(L.REF)
REPORTED = [r["sample_id"] for r in L.samples(reported_only=True)]
BACKGROUND = [r["sample_id"] for r in L.samples() if int(r["sample_no"]) > 10]
LINE = {r["sample_id"]: r["cell_line"] for r in L.samples()}


def tsv(path):
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt") as fh:
        h = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            yield dict(zip(h, line.rstrip("\n").split("\t")))


def bk(s):
    c, p, o = s.rsplit(":", 2)
    return (c, int(p), o)


def a_continuation(A, n):
    """Template orientation: the n reference bases that would follow the A side if the junction moved right
    (i.e. continuing past breakend A away from its retained side)."""
    c, p, o = A
    if o == "+":                                   # template reads ref forward up to p
        return fa.fetch(c, p, p + n)               # ref p+1 .. p+n
    return L.revcomp(fa.fetch(c, max(0, p - 1 - n), p - 1))   # template reads revcomp toward p; continue to p-1..


def b_preceding(B, h, n):
    """Template orientation: the n reference bases that precede the B-side part (which starts h bases inside the
    retained side of B), i.e. what would sit there if the junction moved left."""
    c, p, o = B
    if o == "-":                                   # B part = ref forward from p+h
        return fa.fetch(c, max(0, p + h - 1 - n), p + h - 1)   # ref .. p+h-1
    return L.revcomp(fa.fetch(c, p - h, p - h + n))            # B part = revcomp backward from p-h


junctions = list(tsv(f"{BASE}/stage6_genotype{L.SUFFIX}/junctions.tsv"))
alt = {}
for line in open(f"{BASE}/stage6_genotype{L.SUFFIX}/templates.fa").read().split(">")[1:]:
    name, seq = line.split("\n", 1)
    if name.endswith("|alt"):
        alt[name.split("|")[0]] = seq.replace("\n", "")

# ---- uncertainty (slide) and flank sequences
flank_fa = f"{OUT}/flanks.fa"
ann = {}
with open(flank_fa, "w") as fh:
    for j in junctions:
        jid = j["junction_id"]
        A, B = bk(j["A"]), bk(j["B"])
        T = alt[jid]
        jl, jr = int(j["template_junction_left"]), int(j["template_junction_right"])
        hom, ins = int(j["hom"]), j["ins"]
        x0 = jl + hom                           # template index where the B-side part starts (= len(A part))
        if ins:                                 # an inserted block pins the junction on both sides
            sl = sr = 0
        else:
            bp_ = b_preceding(B, hom, MAXSLIDE)          # T[x0-k:x0] == bp_[-k:]  -> junction can move left by k
            sl = 0
            while sl < MAXSLIDE and x0 - 1 - sl >= 0 and T[x0 - 1 - sl] == bp_[-1 - sl]:
                sl += 1
            ac = a_continuation(A, MAXSLIDE)             # T[x0:x0+k] == ac[:k]    -> junction can move right by k
            sr = 0
            while sr < MAXSLIDE and x0 + sr < len(T) and T[x0 + sr] == ac[sr]:
                sr += 1
        homseq = T[x0 - sl:x0 + sr]
        left_end = x0
        ann[jid] = dict(slide_left=sl, slide_right=sr, uncertainty_bp=sl + sr, microhomology_seq=homseq)
        # flanks exclude the slidable block so that each flank is purely one side of the junction
        fh.write(f">{jid}|A\n{T[max(0, x0 - sl - FL):x0 - sl]}\n>{jid}|B\n{T[jr + sr:jr + sr + FL]}\n")

subprocess.run(f"{ENV}/bwa-mem2 mem -a -t 16 {BWA_IDX} {flank_fa} > {OUT}/flanks_vs_genome.sam 2> {OUT}/flanks_bwa.log",
               shell=True, check=True, executable="/bin/bash")
hits = collections.defaultdict(list)
for r in pysam.AlignmentFile(f"{OUT}/flanks_vs_genome.sam", check_sq=False):
    if r.is_unmapped or r.is_supplementary:
        continue
    hits[r.query_name].append((r.get_tag("AS"), r.reference_name, r.reference_start + 1, r.reference_end,
                               "-" if r.is_reverse else "+", r.mapping_quality, "secondary" if r.is_secondary else "primary"))
with open(f"{OUT}/flank_alignments.tsv", "w") as fh:
    fh.write("flank\tAS\tchrom\tstart\tend\tstrand\tmapq\tclass\n")
    for q, hs in hits.items():
        for h in sorted(hs, reverse=True):
            fh.write(q + "\t" + "\t".join(map(str, h)) + "\n")
for jid in ann:
    for side in ("A", "B"):
        hs = hits.get(f"{jid}|{side}", [])
        bestAS = max((h[0] for h in hs), default=0)
        eq = [h for h in hs if h[0] >= bestAS - MARGIN]
        ann[jid][f"flank{side}_bestAS"] = bestAS
        ann[jid][f"flank{side}_n_equiv"] = len(eq) if hs else 0
        ann[jid][f"flank{side}_competitors"] = ";".join(f"{h[1]}:{h[2]}{h[4]}AS{h[0]}" for h in sorted(eq, reverse=True)[1:6])

# ---- repeat annotation
import bisect
need = collections.defaultdict(list)
for j in junctions:
    for s in ("A", "B"):
        b = bk(j[s])
        need[b[0]].append(b[1])
for c in need:
    need[c].sort()


def load_intervals(fn, cols):
    """keep only intervals that overlap a queried breakend position (UCSC tables are 0-based half-open)."""
    out = collections.defaultdict(list)
    with gzip.open(fn, "rt") as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            c = f[cols[0]]
            if c not in need:
                continue
            s, e = int(f[cols[1]]) + 1, int(f[cols[2]])
            i = bisect.bisect_left(need[c], s)
            if i < len(need[c]) and need[c][i] <= e:
                out[c].append((s, e, cols[3](f)))
    return out


rmsk = segdup = simple = None
if os.path.exists(f"{REFS}/rmsk.txt.gz"):
    # rmsk.txt: 5 genoName, 6 genoStart, 7 genoEnd, 9 repName, 10 repClass, 11 repFamily
    rmsk = load_intervals(f"{REFS}/rmsk.txt.gz", (5, 6, 7, lambda f: f"{f[9]}({f[10]}/{f[11]})"))
    segdup = load_intervals(f"{REFS}/genomicSuperDups.txt.gz", (1, 2, 3, lambda f: f"segdup:{f[7]}:{f[8]}:{float(f[26]):.3f}"))
    simple = load_intervals(f"{REFS}/simpleRepeat.txt.gz", (1, 2, 3, lambda f: f"TRF:{f[16]}x{f[6]}"))


def overlap(db, c, p, pad=0):
    if db is None:
        return "NA"
    hs = [x[2] for x in db.get(c, ()) if x[0] - pad <= p <= x[1] + pad]
    return ";".join(sorted(set(hs))[:5]) or "none"


# ---- support summary
by = collections.defaultdict(dict)
for d in tsv(f"{BASE}/stage6_genotype{L.SUFFIX}/junction_by_sample.tsv"):
    by[d["junction_id"]][d["sample"]] = d
p41 = list(tsv(P41_EVENTS))


def p41_match(A, B):
    for e in p41:
        a = (e["chr1"], int(e["pos1"]))
        b = (e["chr2"], int(e["pos2"]))
        for x, y in ((a, b), (b, a)):
            if x[0] == A[0] and y[0] == B[0] and abs(x[1] - A[1]) <= 15 and abs(y[1] - B[1]) <= 15:
                return e["event_id"], e["max_tier"]
    return "", ""


rows = []
for j in junctions:
    jid = j["junction_id"]
    A, B = bk(j["A"]), bk(j["B"])
    a = ann[jid]
    rep = {s: by[jid][s] for s in REPORTED}
    frag = {s: int(rep[s]["after_frag"]) for s in REPORTED}
    maxf = max(frag.values())
    types = set()
    if any(int(rep[s]["after_split_frag"]) for s in REPORTED):
        types.add("split")
    if any(int(rep[s]["after_disc_frag"]) for s in REPORTED):
        types.add("discordant")
    if int(j["n_contigs"]) > 0:
        types.add("contig")
    bg_n = sum(int(by[jid][s]["after_frag"]) >= 1 for s in BACKGROUND)
    lines = {LINE[s] for s in by[jid] if int(by[jid][s]["after_frag"]) >= 1}
    uniq = a["flankA_n_equiv"] == 1 and a["flankB_n_equiv"] == 1
    if maxf == 0:
        cls = "Not_supported"
    elif maxf == 1:
        cls = "Weak"
    elif not uniq or a["uncertainty_bp"] > 10:
        cls = "Ambiguous_repeat"
    elif maxf >= 3 and len(types) >= 2:
        cls = "Well_supported"
    else:
        cls = "Limited"
    ev, tier = p41_match(A, B)
    in_t = L.in_target(*A[:2]) or L.in_target(*B[:2])
    jtype = ("DEL" if (A[2], B[2]) == ("+", "-") else "DUP" if (A[2], B[2]) == ("-", "+") else "INV") \
        if A[0] == B[0] else "BND"
    rows.append(dict(
        junction_id=jid, in_target=int(in_t), type=jtype, breakend_A=j["A"], breakend_B=j["B"],
        size_bp=(B[1] - A[1]) if A[0] == B[0] else "", classification=cls,
        evidence_types=",".join(sorted(types)), max_after_frag_S01_S10=maxf,
        samples_present_S01_S10=",".join(s.split("_")[0] for s in REPORTED if frag[s] >= 1),
        background_samples_present_S11_S18=bg_n, cell_lines_present=";".join(sorted(lines)),
        microhomology_bp=a["uncertainty_bp"] if not j["ins"] else 0, microhomology_seq=a["microhomology_seq"],
        inserted_bp=j["ins_len"], inserted_seq=j["ins"],
        uncertainty_A=f"{A[1] - a['slide_left'] if A[2] == '+' else A[1] - a['slide_right']}-"
                      f"{A[1] + a['slide_right'] if A[2] == '+' else A[1] + a['slide_left']}",
        uncertainty_B=(lambda h: f"{B[1] + h - a['slide_left']}-{B[1] + h + a['slide_right']}" if B[2] == "-"
                       else f"{B[1] - h - a['slide_right']}-{B[1] - h + a['slide_left']}")(int(j["hom"])),
        uncertainty_bp=a["uncertainty_bp"],
        flankA_n_equiv=a["flankA_n_equiv"], flankB_n_equiv=a["flankB_n_equiv"],
        flankA_competitors=a["flankA_competitors"], flankB_competitors=a["flankB_competitors"],
        repeat_A=overlap(rmsk, A[0], A[1]), repeat_B=overlap(rmsk, B[0], B[1]),
        segdup_A=overlap(segdup, A[0], A[1]), segdup_B=overlap(segdup, B[0], B[1]),
        simple_A=overlap(simple, A[0], A[1]), simple_B=overlap(simple, B[0], B[1]),
        rep_source=j["rep_source"], n_contigs=j["n_contigs"], clusters=j["clusters"],
        P41_event=ev, P41_max_tier=tier,
        **{f"{s.split('_')[0]}_before": rep[s]["before_frag"] for s in REPORTED},
        **{f"{s.split('_')[0]}_after": rep[s]["after_frag"] for s in REPORTED},
        **{f"{s.split('_')[0]}_af": rep[s]["af_W_side"] for s in REPORTED},
        **{f"{s.split('_')[0]}_status": rep[s]["status"] for s in REPORTED}))

with open(f"{OUT}/junction_annotated.tsv", "w") as fh:
    fh.write("\t".join(rows[0].keys()) + "\n")
    for r in rows:
        fh.write("\t".join(map(str, r.values())) + "\n")

# P41 in-target events with and without a P44 counterpart
matched = {r["P41_event"] for r in rows if r["P41_event"]}
with open(f"{OUT}/p41_comparison.tsv", "w") as fh:
    fh.write("P41_event\tP41_max_tier\tP41_in_target\tP44_junction\n")
    m2 = {r["P41_event"]: r["junction_id"] for r in rows if r["P41_event"]}
    for e in p41:
        if e["in_target"] == "True":
            fh.write(f"{e['event_id']}\t{e['max_tier']}\t1\t{m2.get(e['event_id'], '')}\n")
print("stage7 OK", collections.Counter(r["classification"] for r in rows if r["in_target"]))
