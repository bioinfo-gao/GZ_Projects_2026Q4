"""Stage 5 — local reassembly and partner resolution.

(1) Local assembly: for every cluster with >= 2 non-duplicate fragments cohort-wide, reads = all member QNAMEs
    (both mates, full sequences) + reads whose original primary alignment carries a non-adapter clip >= 5 bp
    within TOL of a cluster breakend. Assembled with SPAdes (--only-assembler, k 21,33,55) per cell-line pool
    and per sample (when the sample has >= 2 fragments). Single-fragment clusters are not assembled (a contig
    built from one fragment adds no independent evidence; they are reported as Weak in Stage 7).
(2) Contigs (>= 100 bp) and clip consensus sequences (>= 25 non-N bases) are aligned with `bwa-mem2 mem -a -Y`
    to GRCh38; every alignment (primary, supplementary, secondary) is written to contig_alignments.tsv so that
    competing placements stay visible (plan section 3, repeat-aware localisation).
Outputs: analysis/<RUN>/stage5_assembly{L.SUFFIX}/{asm/<cluster>/<pool>/contigs.fasta, all_contigs.fa, clip_consensus.fa,
         contig_alignments.tsv, contig_junctions.tsv, assembly_log.tsv}
"""
import collections
import gzip
import os
import subprocess
import sys
from concurrent.futures import ThreadPoolExecutor
import pysam
sys.path.insert(0, os.path.dirname(__file__))
import p44lib as L

RUN = os.environ.get("P44_RUN_DATE", "20261024")
BASE = f"{L.PROJ}/analysis/{RUN}_read_rescue"
OUT = f"{BASE}/stage5_assembly{L.SUFFIX}"
os.makedirs(f"{OUT}/asm", exist_ok=True)
ENV = "/Work_bio/gao/configs/.conda/envs/regular_bioinfo/bin"
SPADES = "/Work_bio/gao/configs/.conda/envs/p44_assembly_20261024/bin/spades.py"
BWA_IDX = "/Work_bio/references/Homo_sapiens/GRCh38/GATK.GRCh38/BWAmem2Index/Homo_sapiens_assembly38.fasta"
TOL = 15
SAMPLES = L.samples()
LINE = {r["sample_id"]: r["cell_line"].replace(" ", "") for r in SAMPLES}

clusters = {}
with open(f"{BASE}/stage4_clusters{L.SUFFIX}/clusters.tsv") as fh:
    hdr = fh.readline().rstrip("\n").split("\t")
    for line in fh:
        d = dict(zip(hdr, line.rstrip("\n").split("\t")))
        clusters[d["cluster_id"]] = d
members = collections.defaultdict(list)
with gzip.open(f"{BASE}/stage4_clusters{L.SUFFIX}/cluster_members.tsv.gz", "rt") as fh:
    hdr = fh.readline().rstrip("\n").split("\t")
    for line in fh:
        d = dict(zip(hdr, line.rstrip("\n").split("\t")))
        members[d["cluster_id"]].append(d)


def bk(s):
    if not s:
        return None
    c, p, o = s.rsplit(":", 2)
    return (c, int(p), o)


todo = [c for c, d in clusters.items() if int(d["n_fragments_nondup"]) >= 2]
print("clusters to assemble:", len(todo), "of", len(clusters), flush=True)

# --- per sample: primary read sequences by qname, and local clipped reads near breakends
need = collections.defaultdict(set)          # sample -> qnames
bks_by_sample_need = collections.defaultdict(list)
for c in todo:
    for m in members[c]:
        need[m["sample"]].add(m["qname"])
reads = {}                                     # (sample,qname) -> [(read_end, seq, qual)]
local = collections.defaultdict(lambda: collections.defaultdict(set))   # cluster -> sample -> qnames
bk_list = [(c, bk(clusters[c]["A"])) for c in todo] + [(c, bk(clusters[c]["B"])) for c in todo if clusters[c]["B"]]
for r_s in SAMPLES:
    s = r_s["sample_id"]
    b = pysam.AlignmentFile(f"{BASE}/stage3_rescue/{s}/rescue.bam")
    for r in b.fetch(until_eof=True):
        if r.is_secondary or r.is_supplementary:
            continue
        q = r.query_name
        if not r.is_unmapped and r.cigartuples:
            for side in ("left", "right"):
                op, n = r.cigartuples[0] if side == "left" else r.cigartuples[-1]
                if op != 4 or n < 5 or L.clip_is_adapter(r, side):
                    continue
                pos = r.reference_start + 1 if side == "left" else r.reference_end
                ori = "-" if side == "left" else "+"
                for c, k in bk_list:
                    if k and k[0] == r.reference_name and k[2] == ori and abs(k[1] - pos) <= TOL:
                        local[c][s].add(q)
                        need[s].add(q)
        seq = r.query_sequence
        if seq is None:
            continue
        qual = "".join(chr(x + 33) for x in r.query_qualities) if r.query_qualities is not None else "I" * len(seq)
        if r.is_reverse:
            seq, qual = L.revcomp(seq), qual[::-1]
        reads.setdefault((s, q), {})[1 if r.is_read1 else 2] = (seq, qual)


def write_fq(path, keys):
    n = 0
    with open(path, "w") as fh:
        for k in keys:
            for e, (seq, qual) in sorted(reads.get(k, {}).items()):
                fh.write(f"@{k[0]}|{k[1]}/{e}\n{seq}\n+\n{qual}\n")
                n += 1
    return n


def assemble(job):
    c, pool, keys = job
    d = f"{OUT}/asm/{c}/{pool}"
    os.makedirs(d, exist_ok=True)
    fq = f"{d}/reads.fq"
    n = write_fq(fq, keys)
    if n < 4:
        return (c, pool, len(keys), n, 0, "too_few_reads")
    p = subprocess.run(f"{SPADES} --only-assembler -k 21,33,55 -s {fq} -o {d}/spades -t 2 -m 8 > {d}/spades.log 2>&1",
                       shell=True, executable="/bin/bash")
    cf = f"{d}/spades/contigs.fasta"
    if p.returncode != 0 or not os.path.exists(cf):
        return (c, pool, len(keys), n, 0, f"spades_exit_{p.returncode}")
    out, k = [], 0
    for rec in open(cf).read().split(">")[1:]:
        name, *sq = rec.split("\n")
        sq = "".join(sq)
        if len(sq) >= 100:
            k += 1
            out.append(f">{c}|{pool}|ctg{k}|len{len(sq)}\n{sq}\n")
    with open(f"{d}/contigs.fasta", "w") as fh:
        fh.write("".join(out))
    return (c, pool, len(keys), n, k, "ok")


jobs = []
for c in todo:
    per_s = collections.defaultdict(set)
    for m in members[c]:
        if m["dup"] == "0":
            per_s[m["sample"]].add(m["qname"])
    for s, qs in local[c].items():
        per_s[s] |= qs
    pools = collections.defaultdict(set)
    for s, qs in per_s.items():
        pools[LINE[s]] |= {(s, q) for q in qs}
        if len(qs) >= 2:
            jobs.append((c, s, {(s, q) for q in qs}))
    for line, keys in pools.items():
        jobs.append((c, f"pool_{line}", keys))
print("assembly jobs:", len(jobs), flush=True)
with ThreadPoolExecutor(max_workers=12) as ex:
    logrows = list(ex.map(assemble, jobs))
with open(f"{OUT}/assembly_log.tsv", "w") as fh:
    fh.write("cluster_id\tpool\tn_fragments\tn_reads\tn_contigs_ge100\tstatus\n")
    for r in logrows:
        fh.write("\t".join(map(str, r)) + "\n")

# --- query sequences: contigs + clip consensus
with open(f"{OUT}/all_contigs.fa", "w") as fh:
    for c, pool, *_ in logrows:
        f = f"{OUT}/asm/{c}/{pool}/contigs.fasta"
        if os.path.exists(f):
            fh.write(open(f).read())
with open(f"{OUT}/clip_consensus.fa", "w") as fh:
    for c, d in clusters.items():
        cs = d.get("clip_consensus", "")
        if d["kind"] == "clip_only" and len(cs.replace("N", "")) >= 25:
            a = bk(d["A"])
            fl = 60   # anchor the consensus on 60 bp of its own reference side so placement is relative
            fa = pysam.FastaFile(L.REF)
            if a[2] == "+":
                anchor = fa.fetch(a[0], a[1] - fl, a[1])
                q = anchor + cs
            else:
                anchor = fa.fetch(a[0], a[1] - 1, a[1] - 1 + fl)
                q = cs + anchor
            fh.write(f">{c}|clipcons|anchor{fl}|clip{len(cs)}\n{q}\n")

aln_rows, junc_rows = [], []
for fa_in, tag in ((f"{OUT}/all_contigs.fa", "contig"), (f"{OUT}/clip_consensus.fa", "clipcons")):
    if os.path.getsize(fa_in) == 0:
        continue
    sam = f"{OUT}/{tag}_vs_genome.sam"
    subprocess.run(f"{ENV}/bwa-mem2 mem -a -Y -t 16 {BWA_IDX} {fa_in} > {sam} 2> {OUT}/{tag}_bwa.log",
                   shell=True, check=True, executable="/bin/bash")
    chains = collections.defaultdict(list)
    seqs = {}
    for r in pysam.AlignmentFile(sam, check_sq=False):
        if r.is_unmapped:
            aln_rows.append((tag, r.query_name, "unmapped", "", "", "", "", "", "", "", ""))
            continue
        s = L.seg_from_read(r)
        AS = r.get_tag("AS")
        cls = "secondary" if r.is_secondary else "supplementary" if r.is_supplementary else "primary"
        aln_rows.append((tag, r.query_name, cls, s.chrom, s.start, s.end, s.strand, s.qs, s.qe, r.mapping_quality, AS))
        if cls != "secondary":
            chains[r.query_name].append(s)
        if cls == "primary":
            seqs[r.query_name] = L.revcomp(r.query_sequence) if r.is_reverse else r.query_sequence
    for q, segs in chains.items():
        segs.sort(key=lambda s: s.qs)
        for s1, s2 in zip(segs, segs[1:]):
            j = L.junction_from_segments(s1, s2, seqs.get(q))
            A, B = L.canonical(j["a"], j["b"])
            ins = j["ins"] or ""
            mq1, mq2 = s1.mapq, s2.mapq
            if A != j["a"]:
                ins, mq1, mq2 = L.revcomp(ins), s2.mapq, s1.mapq
            junc_rows.append((tag, q, q.split("|")[0], f"{A[0]}:{A[1]}:{A[2]}", f"{B[0]}:{B[1]}:{B[2]}", j["hom"],
                              j["ins_len"], ins, mq1, mq2))
with open(f"{OUT}/contig_alignments.tsv", "w") as fh:
    fh.write("query_type\tquery\tclass\tchrom\tstart\tend\tstrand\tqstart\tqend\tmapq\tAS\n")
    for r in aln_rows:
        fh.write("\t".join(map(str, r)) + "\n")
with open(f"{OUT}/contig_junctions.tsv", "w") as fh:
    fh.write("query_type\tquery\tcluster_id\tA\tB\thom\tins_len\tins\tmapq_seg1\tmapq_seg2\n")
    for r in junc_rows:
        fh.write("\t".join(map(str, r)) + "\n")
print("stage5 OK", len(aln_rows), "alignments", len(junc_rows), "contig junctions")
