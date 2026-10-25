"""Stage 14b (P45) — WHERE each selected primer binds (P42's 'unique_expected_locus' check, applied to the P45 design).

Stage 14 counts exact reference occurrences; this stage lists them. For every selected primer (primer_design.tsv) and
for the C2 pair that P42 selected (comparison; matched to the P45 event whose two breakends lie within 50 bp of
P42's C2 endpoints), every exact full-length occurrence on either strand of the whole GRCh38 analysis set is located.
  locus_status = junction_specific   0 occurrences (primer crosses the junction; binds only the rearranged allele)
               = expected_locus      1 occurrence, within 1,000 bp of one of the event's breakends
               = off_target          any occurrence elsewhere (reported with coordinates)
Output: analysis/20261025_p45/stage14_primers/{primer_binding_sites.tsv, primer_locus_summary.tsv}
"""
import csv
import os
import sys
import pysam
sys.path.insert(0, os.path.dirname(__file__))
import p45lib as L

OUT = f"{L.OUT}/stage14_primers"
COMP = str.maketrans("ACGTN", "TGCAN")
C2 = (("chr11", 114607191), ("chr21", 20648242))


def rc(s):
    return s.translate(COMP)[::-1]


def tsv(p):
    with open(p) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def bk(s):
    c, p, _ = s.rsplit(":", 2)
    return c, int(p)


design = tsv(f"{OUT}/primer_design.tsv")
events = {r["junction_id"]: r for f in ("junction_table_target.tsv", "junction_table_window_context.tsv")
          for r in tsv(f"{L.OUT}/stage10_tables/{f}")}
queries = []                                  # (source, event_id, side, primer)
for d in design:
    if d["left_primer"]:
        queries += [("P45_selected", d["event_id"], "left", d["left_primer"]),
                    ("P45_selected", d["event_id"], "right", d["right_primer"])]
c2_event = next((j for j, r in events.items()
                 if {bk(r["breakend_A"])[0], bk(r["breakend_B"])[0]} == {"chr11", "chr21"}
                 and all(any(bk(r[k])[0] == c and abs(bk(r[k])[1] - p) <= 50 for k in ("breakend_A", "breakend_B"))
                         for c, p in C2)), "")
for d in tsv(f"{L.P42IN}/results/assay_design_selection.tsv"):
    if d["candidate_label"] == "C2" and d["left_primer"]:
        queries += [("P42_C2_selected", c2_event, "left", d["left_primer"]),
                    ("P42_C2_selected", c2_event, "right", d["right_primer"])]
print("primers:", len(queries), "; C2 matched to event", c2_event or "NONE", flush=True)

pats = {}
for q in queries:
    pats[q[3]] = (q[3].upper(), rc(q[3].upper()))
hits = {p: [] for p in pats}
fa = pysam.FastaFile(L.REF)
for contig in fa.references:
    seq = fa.fetch(contig).upper()
    for p, (fw, rv) in pats.items():
        for strand, s in (("+", fw), ("-", rv)):
            i = seq.find(s)
            while i >= 0:
                hits[p].append((contig, i + 1, strand))
                i = seq.find(s, i + 1)
            if fw == rv:
                break

rows_b, rows_s = [], []
for src, ev, side, p in queries:
    e = events.get(ev)
    ends = [bk(e["breakend_A"]), bk(e["breakend_B"])] if e else []
    for c, pos, st in hits[p]:
        near = any(c == ec and abs(pos - ep) <= 1000 for ec, ep in ends)
        rows_b.append((src, ev, side, p, c, pos, st, int(near)))
    n = len(hits[p])
    n_exp = sum(1 for c, pos, _ in hits[p] if any(c == ec and abs(pos - ep) <= 1000 for ec, ep in ends))
    status = "junction_specific" if n == 0 else ("expected_locus" if n == 1 and n_exp == 1 else "off_target")
    rows_s.append((src, ev, side, p, n, n_exp, status))
with open(f"{OUT}/primer_binding_sites.tsv", "w") as fh:
    fh.write("source\tevent_id\tprimer_side\tprimer\tchrom\tpos_1based\tstrand\twithin_1kb_of_event_breakend\n")
    for r in rows_b:
        fh.write("\t".join(map(str, r)) + "\n")
with open(f"{OUT}/primer_locus_summary.tsv", "w") as fh:
    fh.write("source\tevent_id\tprimer_side\tprimer\texact_reference_sites\tsites_at_event\tlocus_status\n")
    for r in rows_s:
        fh.write("\t".join(map(str, r)) + "\n")
for r in rows_s:
    print("\t".join(map(str, r)))
print("stage14b OK")
