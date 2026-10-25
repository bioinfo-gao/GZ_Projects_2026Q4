"""Stage 4 — de novo junction discovery from the rescue sets of all 18 samples (no P41 candidate gate).

Evidence items (one row each; nothing dropped, every exclusion has a reason):
  split  consecutive aligned segments of one read-end (primary+supplementary chain), in two views:
         'orig' = original BAM records, 'realn' = Stage 3b bwa-mem2 -a re-alignment
  clip   soft clip >= MIN_CLIP bp that is not adapter read-through and not explained by a split segment;
         also taken from realn SECONDARY placements inside W (the explicit form of XA hits)
  disc   read pair (orig primary records) that is not FR with |TLEN| <= MAXINS on one chromosome
Only items with >= 1 breakend inside W are kept.
Counting units (plan section 3): read-end = (sample,qname,read_end); fragment = (sample,qname). Duplicate-flagged
fragments are kept and labelled; counts exclude them.

Clustering (cohort-wide): split items by orientation class and both breakends within TOL bp; discordant pairs attach
to a split cluster when geometrically consistent, otherwise form discordant-only clusters; unexplained clips form
one-sided clusters (partner resolved in Stage 5).
Outputs: analysis/<RUN>/stage4_clusters{L.SUFFIX}/{evidence_items.tsv.gz, clusters.tsv, cluster_members.tsv.gz,
         filter_impact.tsv, maxins.tsv}
"""
import collections
import gzip
import os
import statistics
import sys
import numpy as np
import pysam
sys.path.insert(0, os.path.dirname(__file__))
import p44lib as L

RUN = os.environ.get("P44_RUN_DATE", "20261024")
BASE = f"{L.PROJ}/analysis/{RUN}_read_rescue"
OUT = f"{BASE}/stage4_clusters{L.SUFFIX}"
os.makedirs(OUT, exist_ok=True)
TOL = 15          # bp, split-cluster breakend tolerance (same as P41)
CLIP_TOL = 5      # bp, one-sided clip cluster tolerance
SEC_MARGIN = 10   # AS points; secondary placement must score within this of the primary to count
impact = collections.Counter()


def read_end(r):
    return 1 if r.is_read1 else 2


def orig_seq(r):
    q = r.query_sequence
    return L.revcomp(q) if r.is_reverse else q


def chain_items(sample, recs, view, dupq, zr):
    """split + clip items from one view. recs: list of AlignedSegment."""
    items = []
    by_end = collections.defaultdict(list)
    for r in recs:
        if r.is_unmapped:
            continue
        by_end[(r.query_name, read_end(r))].append(r)
    for (q, e), rs in by_end.items():
        chain = [r for r in rs if not r.is_secondary]
        prim = [r for r in chain if not r.is_supplementary]
        qseq = orig_seq(prim[0]) if prim and prim[0].query_sequence else None
        segs = sorted(((L.seg_from_read(r), r) for r in chain), key=lambda x: x[0].qs)
        covered = [(s.qs, s.qe) for s, _ in segs]
        iend = L.insert_end(qseq, e) if qseq else 10 ** 6
        for (s1, r1), (s2, r2) in zip(segs, segs[1:]):
            # adapter read-through: the 3' segment starts at/after the adapter -> not a genomic junction
            if s2.qs >= iend - 5:
                impact[f"{view}_split_adapter_readthrough"] += 1
                continue
            if qseq and min(L.entropy(qseq[s1.qs:s1.qe]), L.entropy(qseq[s2.qs:min(s2.qe, iend)])) < L.LOW_COMPLEXITY_BITS:
                impact[f"{view}_split_low_complexity_segment"] += 1
                continue
            j = L.junction_from_segments(s1, s2, qseq)
            a, b = j["a"], j["b"]
            if not (L.in_w(*a[:2]) or L.in_w(*b[:2])):
                impact[f"{view}_split_no_end_in_W"] += 1
                continue
            (A, B) = L.canonical(a, b)
            mA, mB = (s1.mapq, s2.mapq) if A == a else (s2.mapq, s1.mapq)
            ins = j["ins"] or ""
            if A != a:                      # junction read from the B side: inserted bases on the other strand
                ins = L.revcomp(ins)
            items.append(dict(sample=sample, qname=q, read_end=e, view=view, type="split", src="chain",
                              A=A, B=B, hom=j["hom"], ins=ins, ins_len=j["ins_len"],
                              mapqA=mA, mapqB=mB, dup=int(q in dupq), zr=zr.get(q, "")))
            impact[f"{view}_split_kept"] += 1
        for r in rs:                                     # clips (chain records and realn secondaries)
            if not L.overlaps_w(r.reference_name, r.reference_start + 1, r.reference_end):
                continue
            # bwa writes secondary records with SEQ='*': take the read-end's primary sequence, re-oriented
            if r.is_secondary:
                # a secondary placement is evidence only if it COMPETES with the primary (AS within MARGIN);
                # partial hits far below the primary score are counted here and not used (plan section 3)
                pAS = prim[0].get_tag("AS") if prim and prim[0].has_tag("AS") else None
                if pAS is None or r.get_tag("AS") < pAS - SEC_MARGIN:
                    impact[f"{view}_secondary_not_competitive"] += 1
                    continue
                impact[f"{view}_secondary_competitive"] += 1
            rseq = r.query_sequence or (None if qseq is None else (L.revcomp(qseq) if r.is_reverse else qseq))
            if rseq is None:
                impact[f"{view}_clip_record_without_sequence"] += 1
                continue
            ops = r.cigartuples
            seg = L.seg_from_read(r)
            if seg.qs >= iend - 5:                       # aligned part lies after the adapter: not genomic
                impact[f"{view}_segment_after_adapter"] += 1
                continue
            for side in ("left", "right"):
                op, n = (ops[0] if side == "left" else ops[-1])
                if op != 4:
                    continue
                if n < L.MIN_CLIP:
                    impact[f"{view}_clip_lt{L.MIN_CLIP}"] += 1
                    continue
                if L.clip_is_adapter(r, side, rseq):
                    impact[f"{view}_clip_adapter"] += 1
                    continue
                # clipped query interval in original orientation
                left_in_read = (side == "left") != r.is_reverse
                tot = len(rseq)
                ci = (0, n) if left_in_read else (tot - n, tot)
                if not left_in_read and iend < tot:          # 3' clip and the read carries adapter
                    if iend <= tot - n + 10:                 # adapter starts within first 10 clip bases
                        impact[f"{view}_clip_adapter_inside"] += 1
                        continue
                    k = iend - (tot - n)                     # genuine insert bases before the adapter
                    if k < L.MIN_CLIP:
                        impact[f"{view}_clip_insert_part_lt{L.MIN_CLIP}"] += 1
                        continue
                    impact[f"{view}_clip_trimmed_to_insert"] += 1
                    n_keep = k
                else:
                    n_keep = n
                if not r.is_secondary:
                    cov = sum(max(0, min(ci[1], qe) - max(ci[0], qs)) for qs, qe in covered)
                    if cov >= 0.8 * n:
                        impact[f"{view}_clip_explained_by_split"] += 1
                        continue
                # clip bases in reference orientation, adjacent to the breakend, trimmed to the insert part
                cseq = rseq[n - n_keep:n] if side == "left" else rseq[len(rseq) - n:len(rseq) - n + n_keep]
                if L.entropy(cseq) < L.LOW_COMPLEXITY_BITS:
                    impact[f"{view}_clip_low_complexity"] += 1
                    continue
                bk = (r.reference_name, seg.start, "-") if side == "left" else (r.reference_name, seg.end, "+")
                if not L.in_w(*bk[:2]):
                    impact[f"{view}_clip_breakend_outside_W"] += 1
                    continue
                cs = cseq
                items.append(dict(sample=sample, qname=r.query_name, read_end=read_end(r), view=view, type="clip",
                                  src="secondary" if r.is_secondary else ("supplementary" if r.is_supplementary
                                                                          else "primary"),
                                  A=bk, B=None, hom=0, ins=cs, ins_len=len(cs), mapqA=r.mapping_quality, mapqB=-1,
                                  dup=int(r.query_name in dupq), zr=zr.get(r.query_name, "")))
                impact[f"{view}_clip_kept_{'secondary' if r.is_secondary else 'chain'}"] += 1
    return items


def disc_items(sample, recs, maxins, dupq, zr):
    items = []
    prim = collections.defaultdict(dict)
    for r in recs:
        if r.is_secondary or r.is_supplementary or r.is_unmapped:
            continue
        prim[r.query_name][read_end(r)] = r
    for q, d in prim.items():
        if len(d) != 2:
            impact["disc_pair_incomplete_or_mate_unmapped"] += 1
            continue
        r1, r2 = d[1], d[2]
        if min(L.entropy(r1.query_sequence or ""), L.entropy(r2.query_sequence or "")) < L.LOW_COMPLEXITY_BITS:
            impact["disc_low_complexity_mate"] += 1
            continue
        s1, s2 = L.seg_from_read(r1), L.seg_from_read(r2)
        if s1.chrom == s2.chrom and s1.strand != s2.strand:
            plus, minus = (s1, s2) if s1.strand == "+" else (s2, s1)
            if plus.start <= minus.end and (minus.end - plus.start + 1) <= maxins:
                impact["pair_concordant"] += 1
                continue
        a = (s1.chrom, s1.end, "+") if s1.strand == "+" else (s1.chrom, s1.start, "-")
        b = (s2.chrom, s2.end, "+") if s2.strand == "+" else (s2.chrom, s2.start, "-")
        if not (L.overlaps_w(s1.chrom, s1.start, s1.end) or L.overlaps_w(s2.chrom, s2.start, s2.end)):
            impact["disc_no_mate_in_W"] += 1
            continue
        A, B = L.canonical(a, b)
        mA, mB = (r1.mapping_quality, r2.mapping_quality) if A == a else (r2.mapping_quality, r1.mapping_quality)
        items.append(dict(sample=sample, qname=q, read_end=0, view="orig", type="disc", src="pair", A=A, B=B,
                          hom=0, ins="", ins_len=0, mapqA=mA, mapqB=mB, dup=int(q in dupq), zr=zr.get(q, "")))
        impact["disc_kept"] += 1
    return items


# ------------------------------------------------------------------ collect
items, maxins_rows = [], []
for row in L.samples():
    s = row["sample_id"]
    D = f"{BASE}/stage3_rescue/{s}"
    orig = list(pysam.AlignmentFile(f"{D}/rescue.bam", check_sq=False).fetch(until_eof=True))
    realn = list(pysam.AlignmentFile(f"{D}/realigned_a.bam", check_sq=False).fetch(until_eof=True))
    dupq = {r.query_name for r in orig if r.is_duplicate}
    zr = {r.query_name: r.get_tag("ZR") for r in orig if r.has_tag("ZR")}
    tl = [abs(r.template_length) for r in orig if r.is_read1 and r.is_proper_pair and not r.is_secondary
          and not r.is_supplementary and r.reference_name == L.CHR and 0 < abs(r.template_length) < 5000]
    maxins = int(max(500, np.percentile(tl, 99.9))) if tl else 500
    maxins_rows.append((s, len(tl), int(np.median(tl)) if tl else "", maxins))
    items += chain_items(s, orig, "orig", dupq, zr)
    items += chain_items(s, realn, "realn", dupq, zr)
    items += disc_items(s, orig, maxins, dupq, zr)
    print(s, "items so far", len(items), "maxins", maxins, flush=True)
MAXINS = max(m[3] for m in maxins_rows)

# ------------------------------------------------------------------ cluster two-sided split items
two = [i for i in items if i["type"] == "split"]
groups = collections.defaultdict(list)
for i in two:
    groups[(i["A"][0], i["A"][2], i["B"][0], i["B"][2])].append(i)
clusters = []
for g, its in groups.items():
    its.sort(key=lambda i: (i["A"][1], i["B"][1]))
    cur = []
    for i in its:
        placed = False
        for c in reversed(cur):              # clusters are created in increasing A order
            if i["A"][1] - c["a"] > TOL:
                break
            if abs(i["A"][1] - c["a"]) <= TOL and abs(i["B"][1] - c["b"]) <= TOL:
                c["m"].append(i)
                placed = True
                break
        if not placed:
            cur.append(dict(a=i["A"][1], b=i["B"][1], m=[i]))
    for c in cur:
        clusters.append(dict(kind="split", g=g, m=c["m"]))


def mode(xs):
    return collections.Counter(xs).most_common(1)[0][0]


for c in clusters:
    c["A"] = (c["g"][0], mode([i["A"][1] for i in c["m"]]), c["g"][1])
    c["B"] = (c["g"][2], mode([i["B"][1] for i in c["m"]]), c["g"][3])

# ------------------------------------------------------------------ attach / cluster discordant pairs


def consistent(bk_disc, bk_clu, slack):
    if bk_disc[0] != bk_clu[0] or bk_disc[2] != bk_clu[2]:
        return False
    if bk_clu[2] == "+":     # read lies left of breakend
        return bk_clu[1] - slack <= bk_disc[1] <= bk_clu[1] + TOL
    return bk_clu[1] - TOL <= bk_disc[1] <= bk_clu[1] + slack


by_group = collections.defaultdict(list)
for c in clusters:
    by_group[c["g"]].append(c)
rest = []
for i in (x for x in items if x["type"] == "disc"):
    hit = None
    for c in by_group.get((i["A"][0], i["A"][2], i["B"][0], i["B"][2]), ()):
        if consistent(i["A"], c["A"], MAXINS) and consistent(i["B"], c["B"], MAXINS):
            hit = c
            break
    if hit:
        hit["m"].append(i)
    else:
        rest.append(i)
dg = collections.defaultdict(list)
for i in rest:
    dg[(i["A"][0], i["A"][2], i["B"][0], i["B"][2])].append(i)
for g, its in dg.items():
    its.sort(key=lambda i: i["A"][1])
    cur = []
    for i in its:
        for c in cur:
            if abs(i["A"][1] - c["m"][-1]["A"][1]) <= MAXINS and abs(i["B"][1] - c["m"][-1]["B"][1]) <= MAXINS:
                c["m"].append(i)
                break
        else:
            cur.append(dict(m=[i]))
    for c in cur:
        ext = (lambda bks, o: max(bks) if o == "+" else min(bks))
        A = (g[0], ext([i["A"][1] for i in c["m"]], g[1]), g[1])
        B = (g[2], ext([i["B"][1] for i in c["m"]], g[3]), g[3])
        clusters.append(dict(kind="disc_only", g=g, m=c["m"], A=A, B=B))

# ------------------------------------------------------------------ clips: attach to two-sided cluster or one-sided cluster
bk_index = collections.defaultdict(list)
for c in clusters:
    if c["kind"] == "split":
        for bk in (c["A"], c["B"]):
            bk_index[(bk[0], bk[2])].append((bk[1], c))
cg = collections.defaultdict(list)
for i in (x for x in items if x["type"] == "clip"):
    hit = None
    for p, c in bk_index.get((i["A"][0], i["A"][2]), ()):
        if abs(p - i["A"][1]) <= TOL:
            hit = c
            break
    if hit:
        hit["m"].append(i)
        impact["clip_attached_to_split_cluster"] += 1
    else:
        cg[(i["A"][0], i["A"][2])].append(i)
for g, its in cg.items():
    its.sort(key=lambda i: i["A"][1])
    cur = []
    for i in its:
        if cur and abs(i["A"][1] - cur[-1]["m"][-1]["A"][1]) <= CLIP_TOL:
            cur[-1]["m"].append(i)
        else:
            cur.append(dict(m=[i]))
    for c in cur:
        A = (g[0], mode([i["A"][1] for i in c["m"]]), g[1])
        clusters.append(dict(kind="clip_only", g=g, m=c["m"], A=A, B=None))


def consensus(clips, ori):
    """clip sequences in reference orientation; '+' clips extend rightwards from the breakend, '-' leftwards."""
    seqs = clips if ori == "+" else [s[::-1] for s in clips]
    n = max(len(s) for s in seqs)
    out = []
    for k in range(n):
        col = [s[k] for s in seqs if len(s) > k]
        if len(col) < max(1, len(seqs) // 2):
            break
        b, cnt = collections.Counter(col).most_common(1)[0]
        out.append(b if cnt / len(col) >= 0.7 else "N")
    s = "".join(out)
    return s if ori == "+" else s[::-1]


# ------------------------------------------------------------------ write
SAMPLES = [r["sample_id"] for r in L.samples()]
fmt = lambda bk: f"{bk[0]}:{bk[1]}:{bk[2]}" if bk else ""
clusters.sort(key=lambda c: (L.chrom_rank(c["A"][0]), c["A"][1]))
with gzip.open(f"{OUT}/cluster_members.tsv.gz", "wt") as fm, open(f"{OUT}/clusters.tsv", "w") as fc:
    fm.write("cluster_id\tsample\tqname\tread_end\tview\ttype\tsrc\tA\tB\thom\tins_len\tins\tmapqA\tmapqB\tdup\trescue_sources\n")
    cols = ["cluster_id", "kind", "A", "B", "in_target", "n_items", "n_fragments_nondup", "n_samples",
            "n_split_fragments", "n_disc_fragments", "n_clip_fragments", "hom_mode", "ins_mode", "clip_consensus",
            "max_mapqA", "max_mapqB"] + [f"{s}_frag" for s in SAMPLES]
    fc.write("\t".join(cols) + "\n")
    for k, c in enumerate(clusters, 1):
        cid = f"P44C{k:05d}"
        c["id"] = cid
        m = c["m"]
        for i in m:
            fm.write("\t".join(map(str, [cid, i["sample"], i["qname"], i["read_end"], i["view"], i["type"], i["src"],
                                         fmt(i["A"]), fmt(i["B"]), i["hom"], i["ins_len"], i["ins"], i["mapqA"],
                                         i["mapqB"], i["dup"], i["zr"]])) + "\n")
        nd = [i for i in m if not i["dup"]]
        frag = lambda t=None: len({(i["sample"], i["qname"]) for i in nd if t is None or i["type"] == t})
        per = collections.defaultdict(set)
        for i in nd:
            per[i["sample"]].add(i["qname"])
        splits = [i for i in m if i["type"] == "split"]
        clips = [i["ins"] for i in m if i["type"] == "clip"]
        intgt = L.in_target(*c["A"][:2]) or (c["B"] is not None and L.in_target(*c["B"][:2]))
        fc.write("\t".join(map(str, [
            cid, c["kind"], fmt(c["A"]), fmt(c["B"]), int(intgt), len(m), frag(), len(per),
            frag("split"), frag("disc"), frag("clip"),
            mode([i["hom"] for i in splits]) if splits else "",
            mode([i["ins"] for i in splits]) if splits else "",
            consensus(clips, c["A"][2]) if (clips and c["kind"] == "clip_only") else "",
            max(i["mapqA"] for i in m), max((i["mapqB"] for i in m), default=-1)] +
            [len(per.get(s, ())) for s in SAMPLES])) + "\n")
with gzip.open(f"{OUT}/evidence_items.tsv.gz", "wt") as fh:
    fh.write("sample\tqname\tread_end\tview\ttype\tsrc\tA\tB\thom\tins_len\tmapqA\tmapqB\tdup\n")
    for i in items:
        fh.write("\t".join(map(str, [i["sample"], i["qname"], i["read_end"], i["view"], i["type"], i["src"],
                                     fmt(i["A"]), fmt(i["B"]), i["hom"], i["ins_len"], i["mapqA"], i["mapqB"],
                                     i["dup"]])) + "\n")
with open(f"{OUT}/filter_impact.tsv", "w") as fh:
    fh.write("item\tcount\n")
    for k in sorted(impact):
        fh.write(f"{k}\t{impact[k]}\n")
with open(f"{OUT}/maxins.tsv", "w") as fh:
    fh.write("sample\tn_proper_pairs_W\tmedian_insert\tmaxins_used\n")
    for r in maxins_rows:
        fh.write("\t".join(map(str, r)) + "\n")
kinds = collections.Counter(c["kind"] for c in clusters)
print("clusters", dict(kinds), "MAXINS", MAXINS)
print("stage4 OK")
