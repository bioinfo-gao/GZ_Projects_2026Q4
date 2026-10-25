"""Stage 1 — audit table of EVERY alignment record stored inside W (chr11:114,599,000-114,630,500), all 18 samples.
No record is filtered; each row carries class labels so later filters can be audited (plan section 3).
Assertions (fail loudly): no secondary records (BAM made without bwa -a), no hard clips (bwa -Y).
Output: analysis/<RUN>/stage1_inventory/<S>.record_audit.tsv.gz and class_summary.tsv
"""
import gzip
import os
import sys
import collections
import pysam
sys.path.insert(0, os.path.dirname(__file__))
import p44lib as L

RUN = os.environ.get("P44_RUN_DATE", "20261024")
OUT = f"{L.PROJ}/analysis/{RUN}_read_rescue/stage1_inventory"
os.makedirs(OUT, exist_ok=True)

COLS = ["sample", "qname", "read_end", "flag", "class", "is_dup", "is_qcfail", "rname", "pos", "end", "mapq",
        "cigar", "seq_len", "softclip_left", "softclip_right", "hardclip", "clip_left_adapter",
        "clip_right_adapter", "rnext", "pnext", "tlen", "proper_pair", "SA", "XA", "NM", "AS", "XS"]


def cls(r):
    if r.is_secondary:
        return "secondary"
    if r.is_supplementary:
        return "supplementary"
    if r.is_unmapped:
        return "unmapped_placed"
    return "primary"


summary = []
for row in L.samples():
    s, bam = row["sample_id"], row["bam_path"]
    b = pysam.AlignmentFile(bam)
    cnt = collections.Counter()
    with gzip.open(f"{OUT}/{s}.record_audit.tsv.gz", "wt") as fh:
        fh.write("\t".join(COLS) + "\n")
        for r in b.fetch(L.CHR, L.W[0] - 1, L.W[1]):
            c = cls(r)
            ops = r.cigartuples or []
            sl = ops[0][1] if ops and ops[0][0] == 4 else 0
            sr = ops[-1][1] if len(ops) > 1 and ops[-1][0] == 4 else 0
            hc = sum(n for op, n in ops if op == 5)
            al = L.clip_is_adapter(r, "left") if sl else False
            ar = L.clip_is_adapter(r, "right") if sr else False
            tag = lambda t: r.get_tag(t) if r.has_tag(t) else ""
            fh.write("\t".join(map(str, [
                s, r.query_name, 1 if r.is_read1 else 2, r.flag, c, int(r.is_duplicate), int(r.is_qcfail),
                r.reference_name, r.reference_start + 1, r.reference_end or "", r.mapping_quality,
                r.cigarstring or "*", len(r.query_sequence or ""), sl, sr, hc, int(al), int(ar),
                r.next_reference_name, r.next_reference_start + 1, r.template_length, int(r.is_proper_pair),
                tag("SA"), tag("XA"), tag("NM"), tag("AS"), tag("XS")])) + "\n")
            cnt[c] += 1
            cnt["dup"] += r.is_duplicate
            cnt["mapq0"] += (r.mapping_quality == 0)
            cnt["hardclip"] += bool(hc)
            cnt["seq_missing"] += (not r.query_sequence)
            cnt["has_SA"] += r.has_tag("SA")
            cnt["has_XA"] += r.has_tag("XA")
            cnt["clip_ge20_nonadapter"] += ((sl >= L.MIN_CLIP and not al) or (sr >= L.MIN_CLIP and not ar))
            cnt["clip_ge20_adapter"] += ((sl >= L.MIN_CLIP and al) or (sr >= L.MIN_CLIP and ar))
            cnt["qnames"] += 0
    summary.append((s, cnt))
    print(s, dict(cnt), flush=True)
    assert cnt["secondary"] == 0, f"{s}: secondary records present — revisit rescue design"
    assert cnt["hardclip"] == 0, f"{s}: hard clips present — SEQ completeness must be audited per record"
    assert cnt["seq_missing"] == 0, f"{s}: records without SEQ"

keys = ["primary", "supplementary", "unmapped_placed", "secondary", "dup", "mapq0", "hardclip", "seq_missing",
        "has_SA", "has_XA", "clip_ge20_nonadapter", "clip_ge20_adapter"]
with open(f"{OUT}/class_summary.tsv", "w") as fh:
    fh.write("sample\t" + "\t".join(keys) + "\n")
    for s, c in summary:
        fh.write(s + "\t" + "\t".join(str(c[k]) for k in keys) + "\n")
print("stage1 OK")
