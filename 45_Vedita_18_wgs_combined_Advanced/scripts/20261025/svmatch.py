"""SV record -> breakend pair conversion and junction match classes, shared by Stage 8x and G2b.
pairs() is copied from P44 8_caller_concordance.py (plus GRIDSS single breakends); classify() is P44's classify()
plus the best two-end distance. P45 fix (2026-10-25): END is parsed from the INFO text (P44's info.get("END")
always fell back to POS because pysam does not expose END in rec.info; rec.stop is unreliable in DELLY genotyped BCFs). TOL (bp) is a module constant set by the caller script before use.
"""
import re

TOL = 100

BND_RE = re.compile(r"([\[\]])([^:\[\]]+):(\d+)([\[\]])")


_END_RE = re.compile(r"(?:^|;)END=(\d+)")


def _info_end(rec):
    """END from the record's INFO text. pysam hides END from rec.info, and rec.stop comes from the BCF rlen field,
    which DELLY writes as 1 for symbolic alleles in genotyped BCFs, so neither can be trusted."""
    m = _END_RE.search(str(rec).split("\t", 8)[7])
    return int(m.group(1)) if m else rec.stop


def pairs(rec):
    """breakend pairs [(A,B)] or ('INS', chrom, pos); P45 addition: ('SGL', chrom, pos) for single breakends"""
    hdr = rec.header.info

    class _I:                                   # header-safe INFO access (pysam raises on undeclared keys)
        def get(self, k, d=None):
            return rec.info[k] if (k in hdr and k in rec.info) else d

        def __contains__(self, k):
            return k in hdr and k in rec.info
    info = _I()
    t = info.get("SVTYPE", "")
    c, p = rec.chrom, rec.pos
    end = info.get("END") or _info_end(rec)   # P45 fix: pysam never exposes END in rec.info (see _info_end)
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
        if alt.startswith(".") or alt.endswith("."):  # P45: GRIDSS single breakend (no partner)
            return [("SGL", c, p)]
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


RANK = {"none": 0, "one_end": 1, "two_end_pos": 2, "two_end": 3}


def classify(J, recpairs):
    """P44 classify(), plus best two-end distance (P45)."""
    A, B = J
    best, dist = "none", None
    for pr in recpairs:
        if pr[0] == "SGL":
            _, c, p = pr
            cls = "one_end" if any(x[0] == c and abs(x[1] - p) <= TOL for x in (A, B)) else "none"
        elif pr[0] == "INS":
            _, c, p = pr
            cls = "two_end" if (A[0] == B[0] == c and abs(A[1] - p) <= TOL and abs(B[1] - p) <= TOL) else \
                  ("one_end" if any(x[0] == c and abs(x[1] - p) <= TOL for x in (A, B)) else "none")
            if cls == "two_end":
                d = max(abs(A[1] - p), abs(B[1] - p))
                dist = d if dist is None else min(dist, d)
        else:
            a, b = pr
            cls = "none"
            for x, y in ((a, b), (b, a)):
                if x[0] == A[0] and y[0] == B[0] and abs(x[1] - A[1]) <= TOL and abs(y[1] - B[1]) <= TOL:
                    cls = "two_end" if (x[2] == A[2] and y[2] == B[2]) else "two_end_pos"
                    if cls == "two_end":
                        d = max(abs(x[1] - A[1]), abs(y[1] - B[1]))
                        dist = d if dist is None else min(dist, d)
                    break
            if cls == "none" and any(z[0] == w[0] and abs(z[1] - w[1]) <= TOL for z in (a, b) for w in (A, B)):
                cls = "one_end"
        if RANK[cls] > RANK[best]:
            best = cls
    return best, dist


