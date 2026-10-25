"""Stage 3 — build the per-sample rescue set R(S) (plan section 3) and fetch every record of every rescued QNAME.

Sources (a QNAME can have several):
  a  a record of the QNAME is stored inside W (incl. unmapped reads placed at a W mate)
  b  a record's mate coordinate (RNEXT/PNEXT) is inside W
  c  an SA segment overlaps W
  d  an XA alternative hit overlaps W
  e  both mates unmapped, and >=1 mate aligns to the W sequence (minimap2 -x sr; >=40 aligned bases,
     identity >=0.95, Shannon entropy of aligned bases >=1.6 bits)
Then all records of these QNAMEs elsewhere in the genome (mates, supplementary segments) are fetched from the
original BAM by position (RNEXT/PNEXT and SA coordinates) with a QNAME filter.

Outputs in analysis/<RUN>/stage3_rescue/<S>/:
  rescue.bam(+bai)          every record of every rescued QNAME; tag ZR:Z:<sources>
  qname_sources.tsv.gz      per QNAME: sources, record counts, completeness checks
  filter_audit.tsv          before/after counts for every step (incl. regex-superset records excluded)
  rescue_R1/R2.fq.gz        primary-record sequences (full reads; no hard clips, Stage 1) for realignment
Usage: 3_collect_rescue_reads.py <sample_id>
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

ENV = "/Work_bio/gao/configs/.conda/envs/regular_bioinfo/bin"
SAMTOOLS, MINIMAP2 = f"{ENV}/samtools", f"{ENV}/minimap2"
RUN = os.environ.get("P44_RUN_DATE", "20261024")
BASE = f"{L.PROJ}/analysis/{RUN}_read_rescue"
S = sys.argv[1]
BAM = {r["sample_id"]: r["bam_path"] for r in L.samples()}[S]
IN = f"{BASE}/stage2_scan/{S}"
OUT = f"{BASE}/stage3_rescue/{S}"
os.makedirs(OUT, exist_ok=True)
audit = []


def run(cmd):
    subprocess.run(cmd, shell=True, check=True, executable="/bin/bash")


def key(r):
    return (r.query_name, r.flag, r.reference_id, r.reference_start, r.cigarstring)


def sources_of(r):
    src = set()
    if r.is_unmapped:
        if L.in_w(r.reference_name, r.reference_start + 1):
            src.add("a")
    elif L.overlaps_w(r.reference_name, r.reference_start + 1, r.reference_end):
        src.add("a")
    if r.next_reference_name == L.CHR and L.in_w(L.CHR, r.next_reference_start + 1):
        src.add("b")
    if r.has_tag("SA") and any(L.overlaps_w(g.chrom, g.start, g.end) for g in L.parse_sa(r.get_tag("SA"))):
        src.add("c")
    if r.has_tag("XA") and any(L.overlaps_w(g.chrom, g.start, g.end) for g in L.parse_xa(r.get_tag("XA"))):
        src.add("d")
    return src


# ---- (a)-(d) from Stage 2 outputs
qsrc = collections.defaultdict(set)
records = {}
n_in = collections.Counter()
for fn in ("window", "scan"):
    with pysam.AlignmentFile(f"{IN}/{fn}.bam") as b:
        header = b.header
        for r in b:
            n_in[fn] += 1
            src = sources_of(r)
            if not src:
                n_in[f"{fn}_regex_superset_or_no_source"] += 1
                continue
            qsrc[r.query_name] |= src
            records[key(r)] = r
for k, v in n_in.items():
    audit.append(("input_records", k, v))
audit.append(("qnames", "after_a-d", len(qsrc)))

# ---- (e) both-unmapped pairs aligned to the W sequence
wfa = f"{OUT}/W.fa"
run(f"{SAMTOOLS} faidx {L.REF} {L.CHR}:{L.W[0]}-{L.W[1]} > {wfa}")
run(f"{SAMTOOLS} fastq -n {IN}/unmapped_pairs.bam 2>/dev/null | gzip > {OUT}/unmapped_pairs.fq.gz")
n_unm = int(subprocess.run(f"{SAMTOOLS} view -c {IN}/unmapped_pairs.bam", shell=True, check=True,
                           capture_output=True, text=True).stdout)
run(f"{MINIMAP2} -ax sr -t 4 --secondary=no {wfa} {OUT}/unmapped_pairs.fq.gz 2>/dev/null "
    f"| {SAMTOOLS} view -F 4 - > {OUT}/unmapped_pairs_vs_W.sam")


def entropy(s):
    c = collections.Counter(s)
    n = len(s)
    return -sum(v / n * math.log2(v / n) for v in c.values()) if n else 0


e_hits = set()
for line in open(f"{OUT}/unmapped_pairs_vs_W.sam"):
    f = line.rstrip("\n").split("\t")
    ops = L.parse_cigar(f[5])
    m = sum(n for n, op in ops if op in "M=X")
    nm = next((int(t[5:]) for t in f[11:] if t.startswith("NM:i:")), 99)
    q = 0
    aligned = []
    for n, op in ops:          # aligned query bases
        if op in "M=X":
            aligned.append(f[9][q:q + n])
        if op in "MIS=X":
            q += n
    if m >= 40 and 1 - nm / max(m, 1) >= 0.95 and entropy("".join(aligned)) >= 1.6:
        e_hits.add(f[0])
audit.append(("unmapped_pairs", "reads_screened", n_unm))
audit.append(("unmapped_pairs", "qnames_hit_W", len(e_hits)))
for q in e_hits:
    qsrc[q].add("e")
with pysam.AlignmentFile(f"{IN}/unmapped_pairs.bam") as b:
    for r in b:
        if r.query_name in e_hits:
            records[key(r)] = r

# ---- fetch remaining records of rescued QNAMEs by position (mates / SA segments outside W)
targets = collections.defaultdict(set)
for r in list(records.values()):
    if r.next_reference_id >= 0 and not L.in_w(r.next_reference_name, r.next_reference_start + 1):
        targets[r.next_reference_name].add(r.next_reference_start + 1)
    if r.has_tag("SA"):
        for g in L.parse_sa(r.get_tag("SA")):
            if not L.overlaps_w(g.chrom, g.start, g.end):
                targets[g.chrom].add(g.start)
bed = f"{OUT}/fetch_positions.bed"
with open(bed, "w") as fh:
    for c in sorted(targets, key=L.chrom_rank):
        for p in sorted(targets[c]):
            fh.write(f"{c}\t{p - 1}\t{p}\n")
qfile = f"{OUT}/rescued_qnames.txt"
with open(qfile, "w") as fh:
    fh.write("\n".join(sorted(qsrc)) + "\n")
audit.append(("fetch", "positions", sum(len(v) for v in targets.values())))
n_f = 0
if targets:
    run(f"{SAMTOOLS} view -@ 2 -b -M -L {bed} -N {qfile} -o {OUT}/fetched.bam {BAM}")
    with pysam.AlignmentFile(f"{OUT}/fetched.bam") as b:
        for r in b:
            n_f += 1
            records.setdefault(key(r), r)
audit.append(("fetch", "records_returned", n_f))

# ---- completeness per QNAME and output
byq = collections.defaultdict(list)
for r in records.values():
    byq[r.query_name].append(r)
rows = []
missing = collections.Counter()
for q, rs in byq.items():
    prim = {1: 0, 2: 0}
    seen = {(r.reference_name, r.reference_start + 1, "-" if r.is_reverse else "+") for r in rs if not r.is_unmapped}
    miss_sa = 0
    for r in rs:
        if not (r.is_secondary or r.is_supplementary):
            prim[1 if r.is_read1 else 2] += 1
        if r.has_tag("SA"):
            for g in L.parse_sa(r.get_tag("SA")):
                if (g.chrom, g.start, g.strand) not in seen:
                    miss_sa += 1
    complete = prim[1] == 1 and prim[2] == 1 and miss_sa == 0
    missing["complete" if complete else "incomplete"] += 1
    rows.append((q, "".join(sorted(qsrc[q])), len(rs), prim[1], prim[2], miss_sa,
                 int(any(r.is_duplicate for r in rs)), int(complete)))
audit.append(("qnames", "total", len(byq)))
audit.append(("qnames", "complete_pairs_and_SA", missing["complete"]))
audit.append(("qnames", "incomplete", missing["incomplete"]))
for s in "abcde":
    audit.append(("qnames_by_source", s, sum(1 for v in qsrc.values() if s in v)))
    audit.append(("qnames_only_source", s, sum(1 for v in qsrc.values() if v == {s})))

with gzip.open(f"{OUT}/qname_sources.tsv.gz", "wt") as fh:
    fh.write("qname\tsources\tn_records\tn_R1_primary\tn_R2_primary\tmissing_SA_segments\tany_dup\tcomplete\n")
    for row in sorted(rows):
        fh.write("\t".join(map(str, row)) + "\n")

tmp = f"{OUT}/rescue.unsorted.bam"
with pysam.AlignmentFile(tmp, "wb", header=header) as o:
    for r in records.values():
        r.set_tag("ZR", "".join(sorted(qsrc[r.query_name])), "Z")
        o.write(r)
run(f"{SAMTOOLS} sort -@ 2 -o {OUT}/rescue.bam {tmp} && {SAMTOOLS} index {OUT}/rescue.bam && rm -f {tmp}")
run(f"{SAMTOOLS} sort -n -@ 2 {OUT}/rescue.bam 2>/dev/null | {SAMTOOLS} fastq -F 0x900 "
    f"-1 {OUT}/rescue_R1.fq.gz -2 {OUT}/rescue_R2.fq.gz -s {OUT}/rescue_single.fq.gz -0 /dev/null - 2>/dev/null")
with open(f"{OUT}/filter_audit.tsv", "w") as fh:
    fh.write("step\titem\tcount\n")
    for a in audit:
        fh.write("\t".join(map(str, a)) + "\n")
open(f"{OUT}/.done", "w").close()
print(S, "\n".join("\t".join(map(str, a)) for a in audit))
