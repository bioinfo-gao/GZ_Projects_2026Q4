"""Stage 6 — final junction set, alt/ref templates, and competitive genotyping in all 18 samples.

Final junctions (precise = both breakends at bp):
  * every Stage 4 split cluster;
  * every Stage 5 contig junction with >= 1 breakend in W (assembly-derived; may refine or add junctions);
  * every clip-consensus junction (one-sided clip cluster whose clipped sequence was placed by bwa-mem2 -a).
  Junctions within TOL bp at both ends (same orientations) are merged; all sources are kept per junction.
  Discordant-only clusters without a contig remain IMPRECISE (counted by pairs only; no bp template).
Templates: alt = A-side 400 bp (read toward the junction) + inserted bases + B-side 400 bp (microhomology kept
once, on the A side); refA / refB = reference +-400 bp around each breakend.
Read-end support (all rescue reads of all samples, bwa-mem2 -a -Y against the template set):
  alt   alignment on alt crosses the junction with >= OVH aligned bases beyond each side of the
        homology/insertion block, AS_alt >= AS_genome + MARGIN (AS_genome = best AS of that read-end in the
        Stage 3b genome re-alignment, primary record), and AS_alt >= AS_ref + MARGIN
  ref   alignment on refA (or refB) crossing the breakend with >= OVH on both sides, AS_ref >= AS_alt and
        AS_ref >= AS_genome - 5 (the read's best genome placement is compatible with this site)
Fragments: alt fragment = QNAME with an alt read-end, or a discordant pair from a linked cluster whose mates sit
on the correct sides within MAXINS of both final breakends. Duplicate-flagged fragments are excluded from counts.
Breakends outside W: ref-spanning read-ends counted directly in the original BAM (rescue reads do not sample them).
'before' = P41-style standard count: original records, non-duplicate split reads with both segments MAPQ >= 20,
plus discordant pairs with both mates MAPQ >= 20, from linked clusters.
Outputs: analysis/<RUN>/stage6_genotype{L.SUFFIX}/{junctions.tsv, templates.fa, junction_by_sample.tsv,
         read_support.tsv.gz, imprecise_disc_clusters.tsv}
"""
import collections
import gzip
import math
import os
import subprocess
import sys
import pysam
sys.path.insert(0, os.path.dirname(__file__))
import p44lib as L

RUN = os.environ.get("P44_RUN_DATE", "20261024")
BASE = f"{L.PROJ}/analysis/{RUN}_read_rescue"
OUT = f"{BASE}/stage6_genotype{L.SUFFIX}"
os.makedirs(OUT, exist_ok=True)
ENV = "/Work_bio/gao/configs/.conda/envs/regular_bioinfo/bin"
FLANK, OVH, MARGIN, TOL = 400, 20, 10, 15
SAMPLES = [r["sample_id"] for r in L.samples()]
BAMS = {r["sample_id"]: r["bam_path"] for r in L.samples()}
fa = pysam.FastaFile(L.REF)


def bk(s):
    if not s:
        return None
    c, p, o = s.rsplit(":", 2)
    return (c, int(p), o)


def fmt(b):
    return f"{b[0]}:{b[1]}:{b[2]}"


def tsv(path):
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt") as fh:
        h = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            yield dict(zip(h, line.rstrip("\n").split("\t")))


clusters = {d["cluster_id"]: d for d in tsv(f"{BASE}/stage4_clusters{L.SUFFIX}/clusters.tsv")}
members = collections.defaultdict(list)
for d in tsv(f"{BASE}/stage4_clusters{L.SUFFIX}/cluster_members.tsv.gz"):
    members[d["cluster_id"]].append(d)
maxins = {d["sample"]: int(d["maxins_used"]) for d in tsv(f"{BASE}/stage4_clusters{L.SUFFIX}/maxins.tsv")}
MAXINS = max(maxins.values())

# ------------------------------------------------------------------ candidate junctions from all sources
cands = []   # dict(A,B,hom,ins,source,ref_id)
for cid, d in clusters.items():
    if d["kind"] == "split":
        cands.append(dict(A=bk(d["A"]), B=bk(d["B"]), hom=int(d["hom_mode"] or 0),
                          ins=d["ins_mode"] if d["ins_mode"] not in ("", "None") else "", source="split", ref=cid))
contig_support = collections.defaultdict(set)
cj = f"{BASE}/stage5_assembly{L.SUFFIX}/contig_junctions.tsv"
if os.path.exists(cj):
    for d in tsv(cj):
        A, B = bk(d["A"]), bk(d["B"])
        if not (L.in_w(*A[:2]) or L.in_w(*B[:2])):
            continue
        src = "contig" if d["query_type"] == "contig" else "clipcons"
        cands.append(dict(A=A, B=B, hom=int(d["hom"]), ins=d["ins"], source=src, ref=d["query"],
                          link=d["cluster_id"]))

# merge
cands.sort(key=lambda c: (c["A"][0], c["A"][2], c["B"][0], c["B"][2], c["A"][1]))
junctions = []
for c in cands:
    for j in junctions[-200:][::-1]:
        if (j["A"][0], j["A"][2], j["B"][0], j["B"][2]) == (c["A"][0], c["A"][2], c["B"][0], c["B"][2]) and \
                abs(j["A"][1] - c["A"][1]) <= TOL and abs(j["B"][1] - c["B"][1]) <= TOL:
            j["sources"].append(c)
            break
    else:
        junctions.append(dict(A=c["A"], B=c["B"], hom=c["hom"], ins=c["ins"], sources=[c]))


def choose_rep(j):
    """representative breakpoint: a contig-derived one if any contig agrees, else the split-cluster mode."""
    for kind in ("contig", "split", "clipcons"):
        s = [c for c in j["sources"] if c["source"] == kind]
        if s:
            r = s[0]
            j.update(A=r["A"], B=r["B"], hom=r["hom"], ins=r["ins"], rep_source=kind)
            return


for k, j in enumerate(junctions, 1):
    choose_rep(j)
    j["id"] = f"P44J{k:04d}"
    j["clusters"] = sorted({c["ref"] for c in j["sources"] if c["source"] == "split"} |
                           {c.get("link") for c in j["sources"] if c.get("link")})
    j["n_contigs"] = len({c["ref"] for c in j["sources"] if c["source"] == "contig"})
    j["clipcons"] = int(any(c["source"] == "clipcons" for c in j["sources"]))

# discordant-only clusters: attach to a precise junction when consistent, else keep as imprecise
def disc_consistent(dA, dB, J, slack):
    def ok(d, b):
        if d[0] != b[0] or d[2] != b[2]:
            return False
        return (b[1] - slack <= d[1] <= b[1] + TOL) if b[2] == "+" else (b[1] - TOL <= d[1] <= b[1] + slack)
    return (ok(dA, J["A"]) and ok(dB, J["B"])) or (ok(dA, J["B"]) and ok(dB, J["A"]))


imprecise = []
for cid, d in clusters.items():
    if d["kind"] != "disc_only":
        continue
    A, B = bk(d["A"]), bk(d["B"])
    hit = [j for j in junctions if disc_consistent(A, B, j, MAXINS)]
    if hit:
        hit[0]["clusters"].append(cid)
    else:
        imprecise.append(d)


# ------------------------------------------------------------------ templates
def side_seq(b, n, toward_junction):
    """sequence of n bp on the retained side of breakend b, oriented so that it reads TOWARD the junction
    (toward_junction=True, used for the A side) or AWAY from it (B side)."""
    c, p, o = b
    if o == "+":
        s = fa.fetch(c, max(0, p - n), p)               # ends at p
        return s if toward_junction else L.revcomp(s)
    s = fa.fetch(c, p - 1, p - 1 + n)                    # starts at p
    return L.revcomp(s) if toward_junction else s


def trim_away(b, h):
    c, p, o = b
    return (c, p + h, o) if o == "-" else (c, p - h, o)


with open(f"{OUT}/templates.fa", "w") as fh:
    for j in junctions:
        left = side_seq(j["A"], FLANK, True)
        right = side_seq(trim_away(j["B"], j["hom"]), FLANK, False)
        j["alt_seq"] = left + j["ins"] + right
        j["jl"] = len(left) - j["hom"]                       # homology block = [jl, jl+hom) on the template
        j["jr"] = len(left) + len(j["ins"])                  # first base after insertion
        fh.write(f">{j['id']}|alt\n{j['alt_seq']}\n")
        for tag, b in (("refA", j["A"]), ("refB", j["B"])):
            c, p, _ = b
            s0 = max(0, p - 1 - FLANK)
            fh.write(f">{j['id']}|{tag}|{s0}\n{fa.fetch(c, s0, p - 1 + FLANK)}\n")
            j[f"{tag}_off"] = p - 1 - s0                     # 0-based index of the breakend base
subprocess.run(f"{ENV}/bwa-mem2 index {OUT}/templates.fa > {OUT}/templates.index.log 2>&1", shell=True, check=True)

# ------------------------------------------------------------------ genome AS per read-end (Stage 3b)
genome_as = {}
for s in SAMPLES:
    for r in pysam.AlignmentFile(f"{BASE}/stage3_rescue/{s}/rescue.bam".replace("rescue.bam", "realigned_a.bam")):
        if r.is_secondary or r.is_supplementary or r.is_unmapped:
            continue
        genome_as[(s, r.query_name, 1 if r.is_read1 else 2)] = r.get_tag("AS")
dupq = collections.defaultdict(set)
for s in SAMPLES:
    for r in pysam.AlignmentFile(f"{BASE}/stage3_rescue/{s}/rescue.bam").fetch(until_eof=True):
        if r.is_duplicate:
            dupq[s].add(r.query_name)

# ------------------------------------------------------------------ align all rescue read-ends to templates
fq = f"{OUT}/all_rescue_reads_SE.fq.gz"
with gzip.open(fq, "wt") as out:
    for s in SAMPLES:
        for e in (1, 2):
            with gzip.open(f"{BASE}/stage3_rescue/{s}/rescue_R{e}.fq.gz", "rt") as f:
                while True:
                    h = f.readline()
                    if not h:
                        break
                    seq, plus, qual = f.readline().rstrip("\n"), f.readline(), f.readline().rstrip("\n")
                    t = L.insert_end(seq, e)             # genotype the genomic insert only (adapter trimmed)
                    if t < 30:
                        continue
                    out.write(f"@{s}|{h[1:].split()[0]}|{e}\n{seq[:t]}\n+\n{qual[:t]}\n")
sam = f"{OUT}/reads_vs_templates.bam"
subprocess.run(f"{ENV}/bwa-mem2 mem -a -Y -t 24 {OUT}/templates.fa {fq} 2> {OUT}/templates.bwa.log "
               f"| {ENV}/samtools view -b -F 4 -o {sam} -", shell=True, check=True, executable="/bin/bash")

J = {j["id"]: j for j in junctions}
best = collections.defaultdict(dict)    # (s,q,e) -> {template_name: (AS, ref_start0, ref_end0)}
for r in pysam.AlignmentFile(sam, check_sq=False):
    key = tuple(r.query_name.split("|"))
    key = (key[0], key[1], int(key[2]))
    t = r.reference_name
    AS = r.get_tag("AS")
    cur = best[key].get(t)
    if cur is None or AS > cur[0]:
        best[key][t] = (AS, r.reference_start, r.reference_end)

support = collections.defaultdict(lambda: collections.defaultdict(lambda: collections.defaultdict(set)))
rows = []
for key, ts in best.items():
    s, q, e = key
    gAS = genome_as.get(key)
    if gAS is None:
        continue
    jids = {t.split("|")[0] for t in ts}
    for jid in jids:
        j = J[jid]
        alt = ts.get(f"{jid}|alt")
        refs = {tag: v for t, v in ts.items() if t.startswith(f"{jid}|ref") for tag in [t.split("|")[1]]}
        best_ref = max((v[0] for v in refs.values()), default=-1)
        if alt and alt[1] <= j["jl"] - OVH and alt[2] >= j["jr"] + OVH and alt[0] >= gAS + MARGIN \
                and alt[0] >= best_ref + MARGIN:
            support[jid][s]["alt"].add(q)
            support[jid][s]["alt_readend"].add((q, e))
            rows.append((jid, s, q, e, "alt", alt[0], gAS, best_ref))
            continue
        for tag, v in refs.items():
            off = j[f"{tag}_off"]
            if v[1] <= off - OVH and v[2] >= off + 1 + OVH and v[0] >= (alt[0] if alt else -1) and v[0] >= gAS - 5:
                support[jid][s][tag].add(q)
                rows.append((jid, s, q, e, tag, v[0], gAS, alt[0] if alt else ""))

# discordant pairs from linked clusters, re-checked against the final breakpoints
before = collections.defaultdict(lambda: collections.defaultdict(set))
for j in junctions:
    for cid in j["clusters"]:
        for m in members.get(cid, ()):
            s, q = m["sample"], m["qname"]
            A, B = bk(m["A"]), bk(m["B"])
            if m["type"] == "disc" and disc_consistent(A, B, j, maxins.get(s, MAXINS)):
                support[j["id"]][s]["alt_disc"].add(q)
                if int(m["mapqA"]) >= 20 and int(m["mapqB"]) >= 20 and m["dup"] == "0":
                    before[j["id"]][s].add(q)
            elif m["type"] == "split" and m["view"] == "orig" and m["dup"] == "0" and \
                    int(m["mapqA"]) >= 20 and int(m["mapqB"]) >= 20:
                before[j["id"]][s].add(q)


def outside_w_ref(s, b):
    """ref-spanning read-ends at a breakend outside W, counted in the original BAM (non-dup, primary)."""
    c, p, _ = b
    n_all = n20 = 0
    with pysam.AlignmentFile(BAMS[s]) as bam:
        for r in bam.fetch(c, p - 1, p):
            if r.is_duplicate or r.is_secondary or r.is_supplementary or r.is_unmapped:
                continue
            if r.reference_start + 1 <= p - OVH and r.reference_end >= p + OVH:
                n_all += 1
                n20 += r.mapping_quality >= 20
                if n_all > 5000:
                    break
    return n_all, n20


def p_detect(depth, af, k=3):
    lam = depth * af
    return 1 - sum(math.exp(-lam) * lam ** i / math.factorial(i) for i in range(k))


# ------------------------------------------------------------------ event groups (alternative representations)
# The same read-ends can support several junctions that differ only in where the partner end is placed (e.g. a
# partner made of a repeat present at many loci). Junctions sharing >= 50% of the smaller alt read-end set are
# one EVENT with alternative partner placements; per-sample counts of the event = union over its members.
alt_ends = collections.defaultdict(set)
alt_as = collections.defaultdict(int)
for jid, s, q, e, cls, AS, gAS, comp in rows:
    if cls == "alt":
        alt_ends[jid].add((s, q, e))
        alt_as[jid] += AS
parent = {j["id"]: j["id"] for j in junctions}


def find(x):
    while parent[x] != x:
        parent[x] = parent[parent[x]]
        x = parent[x]
    return x


by_end = collections.defaultdict(set)
for jid, ends in alt_ends.items():
    for k in ends:
        by_end[k].add(jid)
pairs_seen = set()
for k, js in by_end.items():
    js = sorted(js)
    for a_ in js:
        for b_ in js:
            if a_ < b_ and (a_, b_) not in pairs_seen:
                pairs_seen.add((a_, b_))
                inter = len(alt_ends[a_] & alt_ends[b_])
                if inter >= 0.5 * min(len(alt_ends[a_]), len(alt_ends[b_])):
                    parent[find(a_)] = find(b_)
groups = collections.defaultdict(list)
for j in junctions:
    groups[find(j["id"])].append(j["id"])
rep_of, members_of = {}, {}
for g, ms in groups.items():
    rep = max(ms, key=lambda x: (len(alt_ends[x]), alt_as[x], x))
    members_of[rep] = sorted(ms)
    for m in ms:
        rep_of[m] = rep
for j in junctions:
    j["event_group"] = rep_of[j["id"]]
    j["is_representative"] = int(rep_of[j["id"]] == j["id"])
    j["n_alternative_representations"] = len(members_of[rep_of[j["id"]]]) - 1


def merged(jid, s, key):
    """event-level support for representatives (union over members); own support otherwise."""
    if rep_of[jid] != jid:
        return set(support[jid][s][key])
    out = set()
    for m in members_of[jid]:
        out |= support[m][s][key]
    return out


with open(f"{OUT}/junction_by_sample.tsv", "w") as fh:
    cols = ["junction_id", "sample", "before_frag", "after_frag", "after_split_frag", "after_disc_frag",
            "after_readend", "refA_frag", "refB_frag", "refB_method", "refB_mapq20", "af_W_side", "af_lo95",
            "af_hi95", "local_frag_depth", "p_detect3_af0.10", "p_detect3_af0.25", "p_detect3_af0.50", "status",
            "event_group", "is_representative"]
    fh.write("\t".join(cols) + "\n")
    for j in junctions:
        jid = j["id"]
        for s in SAMPLES:
            sp = collections.defaultdict(set, support[jid][s])
            for key in ("alt", "alt_disc", "alt_readend"):
                sp[key] = merged(jid, s, key)
            bef = set()
            for m in (members_of[jid] if rep_of[jid] == jid else [jid]):
                bef |= before[m][s]
            before[jid][s] = bef
            nd = lambda x: {q for q in x if q not in dupq[s]}
            alt_split, alt_disc = nd(sp["alt"]), nd(sp["alt_disc"])
            alt = alt_split | alt_disc
            rA = nd(sp["refA"]) - alt
            if L.in_w(*j["B"][:2]):
                rB, meth, rB20 = len(nd(sp["refB"]) - alt), "template", ""
            else:
                rB, rB20 = outside_w_ref(s, j["B"])
                meth = "original_bam_spanning_readends"
            wref = len(rA) if L.in_w(*j["A"][:2]) else rB
            n = len(alt) + wref
            af = len(alt) / n if n else float("nan")
            z = 1.96
            if n:
                c0 = (af + z * z / (2 * n)) / (1 + z * z / n)
                hw = z * math.sqrt(af * (1 - af) / n + z * z / (4 * n * n)) / (1 + z * z / n)
                lo, hi = max(0, c0 - hw), min(1, c0 + hw)
            else:
                lo = hi = float("nan")
            status = "present" if len(alt) >= 1 else ("absent" if wref >= 10 else "underpowered")
            fh.write("\t".join(map(str, [
                jid, s, len(nd(before[jid][s])), len(alt), len(alt_split), len(alt_disc), len(sp["alt_readend"]),
                len(rA), rB, meth, rB20, f"{af:.3f}", f"{lo:.3f}", f"{hi:.3f}", n,
                f"{p_detect(n, 0.10):.2f}", f"{p_detect(n, 0.25):.2f}", f"{p_detect(n, 0.50):.2f}", status,
                j["event_group"], j["is_representative"]])) + "\n")

with open(f"{OUT}/junctions.tsv", "w") as fh:
    fh.write("junction_id\tA\tB\thom\tins_len\tins\trep_source\tn_sources\tn_contigs\tclipcons\tclusters\t"
             "template_junction_left\ttemplate_junction_right\tevent_group\tis_representative\t"
             "n_alternative_representations\talternative_representations\n")
    for j in junctions:
        alts = [f"{fmt(J[m]['A'])}~{fmt(J[m]['B'])}" for m in members_of[j["event_group"]] if m != j["id"]]
        fh.write("\t".join(map(str, [j["id"], fmt(j["A"]), fmt(j["B"]), j["hom"], len(j["ins"]), j["ins"],
                                     j["rep_source"], len(j["sources"]), j["n_contigs"], j["clipcons"],
                                     ",".join(j["clusters"]), j["jl"], j["jr"], j["event_group"],
                                     j["is_representative"], j["n_alternative_representations"],
                                     ";".join(alts)])) + "\n")
with gzip.open(f"{OUT}/read_support.tsv.gz", "wt") as fh:
    fh.write("junction_id\tsample\tqname\tread_end\tclass\tAS_template\tAS_genome\tAS_competitor\n")
    for r in rows:
        fh.write("\t".join(map(str, r)) + "\n")
with open(f"{OUT}/imprecise_disc_clusters.tsv", "w") as fh:
    if imprecise:
        fh.write("\t".join(imprecise[0].keys()) + "\n")
        for d in imprecise:
            fh.write("\t".join(d.values()) + "\n")
print("stage6 OK junctions", len(junctions), "imprecise", len(imprecise))
