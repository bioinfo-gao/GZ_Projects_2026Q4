"""Stage 11 — figures for the P44 report (house style: scientific-figure-making SOP, Default_A palette).

Fig1  locus overview: repeat annotation, depth ratio (200 bp bins) with step estimates, target junction events
Fig2  chr11 B-allele frequency, DLD1 (S01-S05) and GM12878 population (S06)
Fig3  read rescue: QNAMEs per sample by rescue source (S01-S10)
Fig4  junction support before vs after rescue, and classification of target events
Fig5  base-level junction view for every target event with >= 2 fragments (one PNG/PDF per event)
Every figure passes the colour assertions; Fig2/Fig3 assert the number of drawn samples.
"""
import collections
import gzip
import os
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pysam
sys.path.insert(0, os.path.dirname(__file__))
sys.path.insert(0, "/home/gao/.claude/skills/scientific-figure-making/assets")
import p44lib as L
from pubfig import (FigureStyle, apply_publication_style, finalize_figure, assert_no_palette_collision,
                    assert_visible_on_white, assert_fill_not_near_black, assert_no_yellow_text, legend_outside,
                    assert_no_text_collisions)

RUN = os.environ.get("P44_RUN_DATE", "20261024")
BASE = f"{L.PROJ}/analysis/{RUN}_read_rescue"
OUT = f"{BASE}/stage11_figures{L.SUFFIX}"
os.makedirs(f"{OUT}/junctions", exist_ok=True)
SUP = f"{BASE}/stage9_supporting"
REP = [r["sample_id"] for r in L.samples(reported_only=True)]
SH = {s: s.split("_")[0] for s in REP}

# ---- Default_A palette (SOP table): areas and strokes are separate roles
SEG = {"target": "#80B1D3", "repeat_array": "#FDB462", "segdup": "#B3DE69"}          # annotation areas
CLASS_FILL = {"Well_supported": "#BC80BD", "Limited": "#FB8072", "Ambiguous_repeat": "#BEBADA",
              "Weak": "#CCEBC5", "Not_supported": "#FCCDE5"}                           # classification areas
STROKE = {"S04": "#8E0152", "S05": "#0000CD", "base": "#8B4500", "line": "#000000", "ci": "#003C30"}
assert_no_palette_collision({"Fig1a annotation": list(SEG.values()), "Fig4 classes": list(CLASS_FILL.values())})
assert_fill_not_near_black(list(SEG.values()) + list(CLASS_FILL.values()))
assert_visible_on_white(list(STROKE.values()), role="line")
assert_no_yellow_text(list(STROKE.values()))
apply_publication_style(FigureStyle(font_size=14, axes_linewidth=1.5))


def tsv(path):
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt") as fh:
        h = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            yield dict(zip(h, line.rstrip("\n").split("\t")))


def bk(s):
    c, p, o = s.rsplit(":", 2)
    return (c, int(p), o)


def w_side(row):
    """the breakend of an event that lies in W (target side)"""
    a, b = bk(row["breakend_A"]), bk(row["breakend_B"])
    return a if L.in_w(*a[:2]) else b


tab = f"{BASE}/stage10_tables{L.SUFFIX}/junction_table_target.tsv"
events = list(tsv(tab))
bounds = {(d["mapq_min"], d["sample"]): d for d in tsv(f"{SUP}/depth_boundary/boundary_estimates.tsv")}
bins = [d for d in tsv(f"{SUP}/depth_boundary/depth_bins_200bp.tsv") if d["mapq_min"] == "20"]

# ================================================================ Fig1 locus overview
X0, X1 = 114_595_000, 114_632_000
rep_intervals = []
with gzip.open(f"{L.PROJ}/refs/simpleRepeat.txt.gz", "rt") as fh:
    for line in fh:
        f = line.split("\t")
        if f[1] == "chr11" and int(f[3]) > X0 and int(f[2]) < X1 and int(f[5]) >= 20:
            rep_intervals.append((int(f[2]) + 1, int(f[3]), f[16].strip(), float(f[6])))
fig, ax = plt.subplots(3, 1, figsize=(14, 10), sharex=True, gridspec_kw={"height_ratios": [1, 2.4, 2.4]})
a0 = ax[0]
a0.add_patch(plt.Rectangle((L.TARGET[0], 0.62), L.TARGET[1] - L.TARGET[0], 0.3, color=SEG["target"], ec="black", lw=0.8))
a0.text(L.TARGET[0], 0.98, "target chr11:114,604,000-114,625,500", va="bottom", fontsize=11)
for s0, s1, unit, cn in rep_intervals:
    a0.add_patch(plt.Rectangle((s0, 0.12), s1 - s0, 0.3, color=SEG["repeat_array"], lw=0))
big = max(rep_intervals, key=lambda x: x[1] - x[0]) if rep_intervals else None
if big:
    a0.text(big[0], -0.05, f"tandem repeat array, {len(big[2])} bp unit x {big[3]:.0f} copies "
            f"({big[0]:,}-{big[1]:,})", va="top", fontsize=11)
a0.set_ylim(-0.45, 1.3)
a0.set_yticks([])
a0.spines["left"].set_visible(False)

a1 = ax[1]
xs = np.array([int(d["bin_start"]) + 100 for d in bins])
base = np.mean([[float(d[s]) for d in bins] for s in ("S01_DLD1_pop", "S02_DLD1_c1_ctrl", "S03_DLD1_c2_ctrl")], axis=0)
for s, col, lab in (("S04_DLD1_c3_dox", STROKE["S04"], "S04 (+Dox clone 3)"),
                    ("S05_DLD1_c4_dox", STROKE["S05"], "S05 (+Dox clone 4)")):
    y = np.array([float(d[s]) for d in bins]) / np.where(base > 0.05, base, np.nan)
    a1.plot(xs, y, color=col, lw=1.6, label=lab)
    b = bounds[("20", s)]
    a1.axvspan(int(b["ci95_lo"]), int(b["ci95_hi"]) + 200, color=col, alpha=0.15, lw=0)
for s in ("S02_DLD1_c1_ctrl", "S03_DLD1_c2_ctrl"):
    y = np.array([float(d[s]) for d in bins]) / np.where(base > 0.05, base, np.nan)
    a1.plot(xs, y, color=STROKE["base"], lw=0.8, alpha=0.7, label="S02, S03 (no-Dox clones)" if s.startswith("S02") else None)
a1.axhline(1, color="black", lw=0.6, ls=":")
a1.axhline(0.5, color="black", lw=0.6, ls=":")
a1.set_ylim(0, 2.2)
a1.set_ylabel("depth / mean(S01-S03)\n(200 bp bins, MAPQ >= 20)")
legend_outside(a1, ncol=3, fontsize=11)

a2 = ax[2]
mk = {"S04": "v", "S05": "^"}
for e in events:
    x = w_side(e)[1]
    y = int(e["max_after_frag_S01_S10"])
    if y == 0:
        continue
    samples = e["samples_present"].split(",") if e["samples_present"] else []
    m = "o"
    if e["DLD1_dox_comparison"].startswith("case_specific"):
        m = mk.get(samples[0], "o") if len(samples) == 1 else "D"
    a2.scatter(x, y, s=60, marker=m, color=CLASS_FILL[e["classification"]], edgecolor="black", lw=0.6, zorder=3)
for cls, c in CLASS_FILL.items():
    if cls != "Not_supported":
        a2.scatter([], [], s=60, color=c, edgecolor="black", lw=0.6, label=cls.replace("_", " "))
a2.scatter([], [], s=60, marker="v", color="white", edgecolor="black", label="S04-only (vs S01-S03)")
a2.scatter([], [], s=60, marker="^", color="white", edgecolor="black", label="S05-only (vs S01-S03)")
a2.set_ylabel("max fragments after rescue\n(any sample S01-S10)")
a2.set_xlabel("chr11 position of the target-side breakend (GRCh38)")
legend_outside(a2, ncol=4, fontsize=10)
a2.set_xlim(X0, X1)
a2.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(lambda v, p: f"{v / 1e6:.3f} Mb"))
finalize_figure(fig, f"{OUT}/P44_Fig1_locus_overview", formats=["png", "pdf"], dpi=300)

# ================================================================ Fig2 chr11 BAF
rows = collections.defaultdict(list)
with gzip.open(f"{SUP}/chr11_baf/baf_sites.tsv.gz", "rt") as fh:
    hdr = {}
    for line in fh:
        f = line.rstrip("\n").split("\t")
        if f[0].startswith("#"):
            hdr[f[0][1:]] = f[2:]
            continue
        rows[f[0]].append(f[1:])
show = [("DLD1-136", "S01_DLD1_pop"), ("DLD1-136", "S02_DLD1_c1_ctrl"), ("DLD1-136", "S03_DLD1_c2_ctrl"),
        ("DLD1-136", "S04_DLD1_c3_dox"), ("DLD1-136", "S05_DLD1_c4_dox"), ("GM12878", "S06_GM12878_pop")]
fig, ax = plt.subplots(len(show), 1, figsize=(14, 2.0 * len(show)), sharex=True)
for k, (line, s) in enumerate(show):
    i = hdr[line].index(s)
    pos = np.array([int(r[0]) for r in rows[line]])
    v = np.array([float(r[1 + i]) if r[1 + i] else np.nan for r in rows[line]])
    ax[k].scatter(pos / 1e6, v, s=0.4, color=STROKE["line"], alpha=0.25, rasterized=True)
    ax[k].axvline(L.TARGET[0] / 1e6, color=STROKE["S04"], lw=1)
    ax[k].set_ylim(-0.02, 1.02)
    ax[k].set_yticks([0, 0.5, 1])
    ax[k].set_ylabel(SH.get(s, s.split("_")[0]), rotation=0, ha="right", va="center")
ax[0].set_title("B-allele fraction at sites heterozygous in the cell line's no-Dox baseline "
                "(DLD1: S01-S03; GM12878: S06); vertical line = target start", fontsize=12, loc="left")
ax[-1].set_xlabel("chr11 position (Mb, GRCh38)")
assert len(ax) == len(show)
finalize_figure(fig, f"{OUT}/P44_Fig2_chr11_BAF", formats=["png", "pdf"], dpi=200)

# ================================================================ Fig3 rescue sources
src = {d["sample"]: d for d in tsv(f"{BASE}/stage10_tables{L.SUFFIX}/rescue_source_summary.tsv")}
cats = [("only in W (a)", lambda d: int(d["qnames_only_a"]), "#80B1D3"),
        ("in W + mate/SA/XA outside", lambda d: int(d["qnames_with_a"]) - int(d["qnames_only_a"]), "#B3DE69"),
        ("mate in W only (b)", lambda d: int(d["qnames_only_b"]), "#FDB462"),
        ("SA segment in W only (c)", lambda d: int(d["qnames_only_c"]), "#FB8072"),
        ("XA alternative hit in W only (d)", lambda d: int(d["qnames_only_d"]), "#BEBADA"),
        ("unmapped pair, aligns to W (e)", lambda d: int(d["qnames_only_e"]), "#FCCDE5")]
assert_fill_not_near_black([c for _, _, c in cats])
fig, axx = plt.subplots(figsize=(13, 6))
ss = [s for s in REP if s in src]
left = np.zeros(len(ss))
for lab, fn, col in cats:
    v = np.array([fn(src[s]) for s in ss])
    axx.barh(range(len(ss)), v, left=left, color=col, edgecolor="black", lw=0.5, label=lab)
    left += v
for i, s in enumerate(ss):
    axx.text(left[i], i, f" {int(left[i]):,}", va="center", fontsize=11)
axx.set_yticks(range(len(ss)))
axx.set_yticklabels([SH[s] for s in ss])
axx.invert_yaxis()
axx.set_xlabel("read pairs in the rescue set")
legend_outside(axx, ncol=2, fontsize=10)
assert len(axx.patches) == len(ss) * len(cats)
finalize_figure(fig, f"{OUT}/P44_Fig3_rescue_sources", formats=["png", "pdf"], dpi=300)

# ================================================================ Fig4 before vs after, classes
fig, axx = plt.subplots(1, 2, figsize=(15, 5.5), gridspec_kw={"width_ratios": [2.2, 1]})
b_ = [sum(int(e[f"{SH[s]}_before"]) >= 1 for e in events) for s in ss]
a_ = [sum(int(e[f"{SH[s]}_after"]) >= 1 for e in events) for s in ss]
a2_ = [sum(int(e[f"{SH[s]}_after"]) >= 2 for e in events) for s in ss]
x = np.arange(len(ss))
for off, v, col, lab in ((-0.27, b_, "#80B1D3", "before rescue: >= 1 fragment (standard count, MAPQ >= 20)"),
                         (0.0, a_, "#FDB462", "after rescue: >= 1 fragment"),
                         (0.27, a2_, "#B3DE69", "after rescue: >= 2 fragments")):
    bars = axx[0].bar(x + off, v, width=0.26, color=col, edgecolor="black", lw=0.5, label=lab)
    for bb, vv in zip(bars, v):
        axx[0].text(bb.get_x() + bb.get_width() / 2, vv, str(vv), ha="center", va="bottom", fontsize=9)
axx[0].set_xticks(x)
axx[0].set_xticklabels([SH[s] for s in ss])
axx[0].set_ylabel("target-region junction events")
legend_outside(axx[0], ncol=1, fontsize=10)
cc = collections.Counter(e["classification"] for e in events)
order = ["Well_supported", "Limited", "Ambiguous_repeat", "Weak", "Not_supported"]
vals = [cc.get(k, 0) for k in order]
bars = axx[1].barh(range(len(order)), vals, color=[CLASS_FILL[k] for k in order], edgecolor="black", lw=0.5)
for i, v in enumerate(vals):
    axx[1].text(v, i, f" {v}", va="center", fontsize=11)
axx[1].set_yticks(range(len(order)))
axx[1].set_yticklabels([k.replace("_", " ") for k in order])
axx[1].invert_yaxis()
axx[1].set_xlabel("target-region junction events")
axx[1].set_xlim(0, max(vals) * 1.25 + 1)
finalize_figure(fig, f"{OUT}/P44_Fig4_before_after_and_classes", formats=["png", "pdf"], dpi=300)

# ================================================================ Fig5 per-event base-level junction view
alt_seq = {}
for rec in open(f"{BASE}/stage6_genotype{L.SUFFIX}/templates.fa").read().split(">")[1:]:
    name, seq = rec.split("\n", 1)
    if name.endswith("|alt"):
        alt_seq[name.split("|")[0]] = seq.replace("\n", "")
jinfo = {d["junction_id"]: d for d in tsv(f"{BASE}/stage6_genotype{L.SUFFIX}/junctions.tsv")}
support = collections.defaultdict(list)
for d in tsv(f"{BASE}/stage6_genotype{L.SUFFIX}/read_support.tsv.gz"):
    if d["class"] == "alt":
        support[d["junction_id"]].append(d)
tpl_bam = f"{BASE}/stage6_genotype{L.SUFFIX}/reads_vs_templates.bam"
aln = collections.defaultdict(list)
want = {e["junction_id"] for e in events if int(e["max_after_frag_S01_S10"]) >= 2}
members = collections.defaultdict(set)
for jid in want:
    for m in [jid] + [x for x in jinfo if jinfo[x]["event_group"] == jid and x != jid]:
        members[jid].add(m)
keys = {(d["sample"], d["qname"], d["read_end"]) for jid in want for d in support[jid]}
readseq = {}                                     # original-orientation sequences (secondaries carry SEQ='*')
with gzip.open(f"{BASE}/stage6_genotype{L.SUFFIX}/all_rescue_reads_SE.fq.gz", "rt") as fh:
    while True:
        h = fh.readline()
        if not h:
            break
        sq = fh.readline().strip()
        fh.readline()
        fh.readline()
        k = tuple(h[1:].strip().split("|"))
        if k in keys:
            readseq[k] = sq
for r in pysam.AlignmentFile(tpl_bam, check_sq=False):
    t = r.reference_name.split("|")
    if t[1] != "alt" or t[0] not in want:
        continue
    s, q, e = r.query_name.split("|")
    if (s, q, e) in keys:
        sq = readseq.get((s, q, e))
        if sq is None:
            continue
        aln[t[0]].append((s, q, e, r.reference_start, None, r.get_tag("AS"),
                          r.cigartuples, L.revcomp(sq) if r.is_reverse else sq))
VIEW = 70
fig_paths = []
for e in events:
    jid = e["junction_id"]
    if jid not in want:
        continue
    T = alt_seq[jid]
    jl, jr = int(jinfo[jid]["template_junction_left"]), int(jinfo[jid]["template_junction_right"])
    hom = int(jinfo[jid]["hom"])
    x0 = jl + hom
    v0, v1 = x0 - VIEW, x0 + VIEW
    ok_ends = {(d["sample"], d["qname"], d["read_end"]) for d in support[jid]}
    best = {}
    for s, q, en, st, qa, AS, cig, qs in aln[jid]:
        if (s, q, en) in ok_ends and ((s, q, en) not in best or AS > best[(s, q, en)][2]):
            best[(s, q, en)] = (st, cig, AS, qs)
    readrows = []
    for (s, q, en), (st, cig, AS, qs) in sorted(best.items()):
        line = [" "] * (v1 - v0)
        rp, qp = st, 0
        for op, n in cig:
            if op in (0, 7, 8):
                for k in range(n):
                    if v0 <= rp + k < v1:
                        b = qs[qp + k]
                        line[rp + k - v0] = "." if b == T[rp + k] else b
                rp += n
                qp += n
            elif op in (1, 4):
                qp += n
            elif op in (2, 3):
                for k in range(n):
                    if v0 <= rp + k < v1:
                        line[rp + k - v0] = "-"
                rp += n
        readrows.append((s, "".join(line)))
    n = len(readrows)
    fig, axx = plt.subplots(figsize=(15, 1.6 + 0.26 * (n + 3)))
    axx.set_xlim(-14, v1 - v0 + 1)
    axx.set_ylim(-(n + 2.5), 2.6)
    axx.axis("off")
    ref = T[v0:v1]
    for k, ch in enumerate(ref):
        colr = STROKE["S04"] if k < x0 - hom - v0 else (STROKE["ci"] if k < jr - v0 else STROKE["S05"])
        axx.text(k + 0.5, 0, ch, family="monospace", fontsize=8, ha="center", va="center", color=colr)
    axx.text(-1, 0, "junction", ha="right", va="center", fontsize=9)
    axx.axvline(x0 - v0, color="black", lw=0.6, ls="--")
    A_, B_ = e["breakend_A"], e["breakend_B"]
    axx.text(0, 1.6, f"{jid}  {e['type']}  A={A_}  B={B_}  class={e['classification'].replace('_', ' ')}  "
             f"microhomology={e['microhomology_seq'] or '-'}  inserted={e['inserted_seq'] or '-'}",
             fontsize=9, va="center")
    alt_n = e["n_alternative_representations"]
    partner = (f"B-side sequence class: {e['seqclass_flank_B']}; {alt_n} alternative B placements"
               if alt_n != "0" else f"B-side sequence class: {e['seqclass_flank_B']}; B placement unique")
    axx.text(0, 0.9, "wine = A side, dark teal = microhomology/inserted bases, blue = B side; dashed line = junction; "
             "'.' = read base identical to the junction sequence; " + partner, fontsize=8, va="center")
    for i, (s, row) in enumerate(readrows, 1):
        axx.text(-1, -i, SH.get(s, s.split("_")[0]), ha="right", va="center", fontsize=8)
        for k, ch in enumerate(row):
            if ch != " ":
                axx.text(k + 0.5, -i, ch, family="monospace", fontsize=8, ha="center", va="center",
                         color="black" if ch in ".-" else STROKE["S04"])
    p = finalize_figure(fig, f"{OUT}/junctions/P44_Fig5_{jid}", formats=["png", "pdf"], dpi=200)
    fig_paths.append((jid, n))
with open(f"{OUT}/junctions/index.tsv", "w") as fh:
    fh.write("junction_id\tsupporting_read_ends_drawn\n")
    for jid, n in fig_paths:
        fh.write(f"{jid}\t{n}\n")
print("stage11 OK; per-event figures:", len(fig_paths))
