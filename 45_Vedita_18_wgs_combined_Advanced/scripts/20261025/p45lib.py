"""P45 shared helpers (copied from P44 p44lib.py; only PROJ, SAMPLE_TABLE, W_PAD default and the client-ID
helpers below differ): sample table, locus constants, alignment geometry, SA/XA parsing, adapter test.

Coordinate conventions (used in every stage):
  * positions are 1-based, inclusive (VCF/SAM style).
  * a breakend is (chrom, pos, ori). ori '+' = retained sequence lies to the LEFT of pos (ends at pos, the
    junction is immediately right of pos); ori '-' = retained sequence lies to the RIGHT (starts at pos).
  * query coordinates are 0-based half-open in the ORIGINAL read orientation (as sequenced).
"""
import csv
import os
import re
from dataclasses import dataclass

PROJ = "/home/gao/projects_2026Q4/45_Vedita_18_wgs_combined_Advanced"
SAMPLE_TABLE = f"{PROJ}/scripts/0_sample_table.csv"
REF = "/Work_bio/references/Homo_sapiens/GRCh38/GATK.GRCh38/Homo_sapiens_assembly38.fasta"
CHR = "chr11"
TARGET = (114604000, 114625500)
W_PAD = 50000            # = P44 round 2 (user D3, 2026-10-25); P45 inherits those stage 1-8 products
W = (TARGET[0] - W_PAD, TARGET[1] + W_PAD)
READ_LEN = 150
MIN_CLIP = 20             # soft clip length to count as breakpoint evidence (plan section 3)
ADAPTERS = {               # Illumina TruSeq read-through adapters (identified from data, Stage 1 note)
    "R1": "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA",
    "R2": "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT",
}

COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(s):
    return s.translate(COMP)[::-1]


import os
# Test mode only (defaults = canonical run): P44_SAMPLES=comma list restricts samples; P44_STAGE_SUFFIX
# (e.g. "_test") redirects stage4-10 directories so test output never mixes with canonical output.
SUFFIX = ""             # P45 has no test mode
_ONLY = []


def samples(reported_only=False):
    with open(SAMPLE_TABLE) as fh:
        rows = list(csv.DictReader(fh))
    if _ONLY:
        rows = [r for r in rows if r["sample_id"] in _ONLY]
    if reported_only:
        rows = [r for r in rows if int(r["sample_no"]) <= 10]
    return rows


def in_w(chrom, pos):
    return chrom == CHR and W[0] <= pos <= W[1]


def overlaps_w(chrom, start, end):
    return chrom == CHR and start <= W[1] and end >= W[0]


def in_target(chrom, pos, pad=50):
    return chrom == CHR and TARGET[0] - pad <= pos <= TARGET[1] + pad


# ---------------------------------------------------------------- CIGAR geometry
CIG_RE = re.compile(r"(\d+)([MIDNSHP=X])")


def parse_cigar(cig):
    return [(int(n), op) for n, op in CIG_RE.findall(cig)]


def cigar_ref_len(ops):
    return sum(n for n, op in ops if op in "MDN=X")


def cigar_query_span(ops):
    """(clip_left, aligned_query_len, clip_right, total_len) in CIGAR (reference) orientation; H counted."""
    left = 0
    i = 0
    while i < len(ops) and ops[i][1] in "SH":
        left += ops[i][0]
        i += 1
    right = 0
    j = len(ops) - 1
    while j >= 0 and ops[j][1] in "SH":
        right += ops[j][0]
        j -= 1
    aligned = sum(n for n, op in ops[i:j + 1] if op in "MI=X")
    return left, aligned, right, left + aligned + right


@dataclass
class Seg:
    chrom: str
    start: int      # 1-based first aligned ref base
    end: int        # 1-based last aligned ref base
    strand: str     # '+' / '-'
    qs: int         # query start, original read orientation
    qe: int         # query end (exclusive)
    mapq: int
    src: str        # 'primary' / 'supplementary' / 'secondary' / 'SA' / 'XA'
    nm: int = -1


def seg_from_cigar(chrom, pos, strand, cig, mapq, src, nm=-1):
    ops = parse_cigar(cig)
    cl, al, cr, tot = cigar_query_span(ops)
    end = pos + cigar_ref_len(ops) - 1
    if strand == "+":
        qs, qe = cl, cl + al
    else:
        qs, qe = cr, cr + al
    return Seg(chrom, pos, end, strand, qs, qe, mapq, src, nm)


def seg_from_read(r):
    src = "secondary" if r.is_secondary else "supplementary" if r.is_supplementary else "primary"
    nm = r.get_tag("NM") if r.has_tag("NM") else -1
    return seg_from_cigar(r.reference_name, r.reference_start + 1, "-" if r.is_reverse else "+",
                          r.cigarstring, r.mapping_quality, src, nm)


def parse_sa(tag):
    """SA:Z:chr,pos,strand,CIGAR,mapQ,NM;..."""
    out = []
    for item in tag.strip(";").split(";"):
        if not item:
            continue
        c, p, s, cig, mq, nm = item.split(",")
        out.append(seg_from_cigar(c, int(p), s, cig, int(mq), "SA", int(nm)))
    return out


def parse_xa(tag):
    """XA:Z:chr,+pos,CIGAR,NM;...  (bwa alternative hits; pos signed by strand)"""
    out = []
    for item in tag.strip(";").split(";"):
        if not item:
            continue
        c, sp, cig, nm = item.split(",")
        strand = "-" if sp[0] == "-" else "+"
        out.append(seg_from_cigar(c, int(sp.lstrip("+-")), strand, cig, 0, "XA", int(nm)))
    return out


def junction_from_segments(s1, s2, qseq):
    """Junction between consecutive query segments s1 (earlier in read) and s2.
    Returns dict with breakend A (from s1), breakend B (from s2), microhomology length (query overlap),
    inserted sequence (query gap).  qseq is the read in ORIGINAL orientation (may be None)."""
    a = (s1.chrom, s1.end, "+") if s1.strand == "+" else (s1.chrom, s1.start, "-")
    b = (s2.chrom, s2.start, "-") if s2.strand == "+" else (s2.chrom, s2.end, "+")
    overlap = s1.qe - s2.qs
    hom = max(overlap, 0)
    ins = qseq[s1.qe:s2.qs] if (qseq and overlap < 0) else ("" if overlap >= 0 else None)
    return dict(a=a, b=b, hom=hom, ins=ins, ins_len=max(-overlap, 0))


def canonical(a, b):
    """Order two breakends so that the junction key is strand-independent."""
    ka = (chrom_rank(a[0]), a[1])
    kb = (chrom_rank(b[0]), b[1])
    return (a, b) if ka <= kb else (b, a)


_MAIN = {f"chr{i}": i for i in range(1, 23)}
_MAIN.update({"chrX": 23, "chrY": 24, "chrM": 25})


def chrom_rank(c):
    return (_MAIN.get(c, 100), c)


# ---------------------------------------------------------------- adapter test
def _prefix_match(s, adapter, max_offset=3, max_mm=1, k=13):
    for off in range(0, max_offset + 1):
        a = adapter[off:off + k]
        t = s[:len(a)]
        if len(t) < 6:
            return False
        mm = sum(1 for x, y in zip(t, a) if x != y)
        if mm <= max_mm:
            return True
    return False


def clip_is_adapter(r, side, seq=None):
    """side = 'left'/'right' clip in reference orientation. A 3'-end clip (right on + strand, left on -)
    is adapter read-through if its sequence starts like the TruSeq adapter, or its length equals
    READ_LEN - |TLEN| (+-3) for a short-insert proper pair."""
    three_prime = (side == "right") != r.is_reverse
    if not three_prime:
        return False
    ops = r.cigartuples
    n = ops[-1][1] if side == "right" else ops[0][1]
    q = seq if seq is not None else r.query_sequence
    if q is None:
        return False
    if side == "right":
        clip = q[-n:]
    else:
        clip = revcomp(q[:n])          # bring to read orientation, 5'->3'
    ad = ADAPTERS["R1" if r.is_read1 else "R2"]
    if _prefix_match(clip, ad):
        return True
    tl = abs(r.template_length)
    if r.is_proper_pair and 0 < tl < READ_LEN and abs(n - (READ_LEN - tl)) <= 3:
        return True
    return False


def entropy(seq):
    """Shannon entropy (bits) of base composition; poly-G tails ~0, random DNA ~2."""
    import math
    seq = seq.upper()
    n = len(seq)
    if not n:
        return 0.0
    c = {b: seq.count(b) for b in set(seq)}
    return -sum(v / n * math.log2(v / n) for v in c.values())


LOW_COMPLEXITY_BITS = 1.2


def insert_end(qseq, read_end):
    """0-based start of adapter read-through in a read (ORIGINAL orientation), or len(qseq) if none.
    Scans for the first 13 adapter bases with <= 1 mismatch (also allows the adapter to start 1-3 bases
    in, matching the clipped-adapter rule)."""
    if not qseq:
        return 0
    ad = ADAPTERS["R1" if read_end == 1 else "R2"][:13]
    for i in range(0, len(qseq) - 8):
        t = qseq[i:i + 13]
        if sum(1 for x, y in zip(t, ad) if x != y) + (13 - len(t)) <= 1:
            return i
    return len(qseq)


# ---------------------------------------------------------------- P45 paths and client-facing identifiers
INH = f"{PROJ}/analysis/20261025_inherited"       # copied, read-only products (Stage 0 / 0b)
P44R2 = f"{INH}/P44_r2"                             # P44 round-2 stage1-8 (W +/-50 kb)
SUP = f"{INH}/P44_r1_stage9_supporting"             # depth boundary, chr11 BAF/SNV, VCF package (W-independent)
P42IN = f"{INH}/P42"
OUT = f"{PROJ}/analysis/20261025_p45"


def cid(x):
    """Client-facing junction ID: the internal 'P44J0332' becomes 'J0332' (no internal project numbers in deliverables)."""
    return re.sub(r"\bP44J(\d+)", r"J\1", x) if x else x
