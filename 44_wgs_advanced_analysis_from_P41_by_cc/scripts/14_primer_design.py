"""Stage 14 (P4) — in-silico junction-spanning PCR primers for every event with >= 2 fragments (S01-S10).

Design: Primer3 (primer3-py) on the event's junction template (Stage 6, 400 bp + junction + 400 bp); the target
interval covers the microhomology/inserted block + 5 bp each side, so every product spans the junction.
Product 150-350 bp; primer length 18-27, Tm 57-63 C; up to 50 candidate pairs ranked by Primer3 penalty, in two passes:
"flank" (primers anywhere outside the target) then "junction_overlap" (one primer must cross the junction with >= 4 bases
at its 3' side, so its 3' end exists only in the rearranged allele). Pairs whose two primers are identical or reverse
complements or share any 12-mer (palindromic fold-back template) are discarded. In the junction_overlap pass the
TARGET constraint is dropped (Primer3 never places a primer inside a TARGET); the overlap requirement ensures the product
spans the junction.
Specificity (whole GRCh38 analysis set, both strands, exact string search, one grep -F per query):
  full_hits  = occurrences of the full primer;  end15_hits = occurrences of its 3'-terminal 15 bases.
Strict pass: both primers full_hits <= 1 (0 = the primer crosses the junction, i.e. junction-specific sequence absent from
the reference) and end15_hits <= 1 (3' end primes at no more than one reference site).
Junction-specific pass (second tier): one primer has full_hits 0 and end15_hits <= 1 (cannot prime the reference),
the other exactly one full-length reference site; reported with all counts. Strict pairs are preferred. The first pair (lowest penalty) in which both primers pass
is selected. If none passes the event is reported as "no_unique_pair" with the reason (e.g. partner side is
telomeric repeat, or flank lies in the 42 bp tandem array), never forced.
The genome-side hit of a correct junction primer is its own flank locus (A or B side), so 1 is the expected count.
Output: analysis/<RUN>/stage14_primers/{primer_design.tsv, primer_candidates.tsv}
Run with the P44 env python (primer3-py): /Work_bio/gao/configs/.conda/envs/p44_assembly_20261024/bin/python
"""
import csv
import gzip
import os
import sys
import primer3
import collections
sys.path.insert(0, os.path.dirname(__file__))
import p44lib as L

RUN = os.environ.get("P44_RUN_DATE", "20261024")
BASE = f"{L.PROJ}/analysis/{RUN}_read_rescue"
OUT = f"{BASE}/stage14_primers{L.SUFFIX}"
os.makedirs(OUT, exist_ok=True)
COMP = str.maketrans("ACGTN", "TGCAN")


def rc(s):
    return s.translate(COMP)[::-1]


def tsv(p):
    with open(p) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


events = [r for f in ("junction_table_target.tsv", "junction_table_window_context.tsv")
          for r in tsv(f"{BASE}/stage10_tables{L.SUFFIX}/{f}") if int(r["max_after_frag_S01_S10"]) >= 2]
jinfo = {d["junction_id"]: d for d in tsv(f"{BASE}/stage6_genotype{L.SUFFIX}/junctions.tsv")}
alt = {}
for rec in open(f"{BASE}/stage6_genotype{L.SUFFIX}/templates.fa").read().split(">")[1:]:
    name, seq = rec.split("\n", 1)
    if name.endswith("|alt"):
        alt[name.split("|")[0]] = seq.replace("\n", "").upper()

def related(a, b, k=12):
    """True if the two primers share any k-mer in either orientation (identical / shifted / reverse-complement)."""
    km = {a[i:i + k] for i in range(len(a) - k + 1)}
    return any(b[i:i + k] in km or rc(b)[i:i + k] in km for i in range(len(b) - k + 1))


cands, explain = {}, {}
n_palin = collections.Counter()
for e in events:
    jid = e["junction_id"]
    T = alt[jid]
    ji = jinfo[jid]
    x0 = int(ji["template_junction_left"])                  # start of microhomology block
    x1 = int(ji["template_junction_right"]) + int(ji["hom"])  # end of homology/insert block (exclusive, approx.)
    tgt = (max(0, x0 - 5), (x1 - x0) + 10)
    junction_pos = int(ji["template_junction_right"])          # first base after the junction block
    cands[jid] = []
    for mode in ("flank", "junction_overlap"):
        seq_args = {"SEQUENCE_ID": jid, "SEQUENCE_TEMPLATE": T, "SEQUENCE_TARGET": list(tgt)}
        glob_args = {"PRIMER_NUM_RETURN": 50, "PRIMER_PRODUCT_SIZE_RANGE": [[150, 350]], "PRIMER_OPT_SIZE": 22,
                     "PRIMER_MIN_SIZE": 18, "PRIMER_MAX_SIZE": 27, "PRIMER_OPT_TM": 60.0, "PRIMER_MIN_TM": 57.0,
                     "PRIMER_MAX_TM": 63.0, "PRIMER_MIN_GC": 30.0, "PRIMER_MAX_GC": 70.0, "PRIMER_MAX_POLY_X": 4}
        if mode == "junction_overlap":         # one primer must cross the junction (>= 4 bases at its 3' side);
            del seq_args["SEQUENCE_TARGET"]    # Primer3 never places a primer inside a TARGET, so drop it here
            seq_args["SEQUENCE_OVERLAP_JUNCTION_LIST"] = [junction_pos]
            glob_args["PRIMER_MIN_3_PRIME_OVERLAP_OF_JUNCTION"] = 4
            glob_args["PRIMER_MIN_5_PRIME_OVERLAP_OF_JUNCTION"] = 7
        res = primer3.bindings.design_primers(seq_args, glob_args)
        n = res.get("PRIMER_PAIR_NUM_RETURNED", 0)
        for i in range(n):
            l, r = res[f"PRIMER_LEFT_{i}_SEQUENCE"], res[f"PRIMER_RIGHT_{i}_SEQUENCE"]
            if related(l, r):                   # palindromic template (short fold-back): single-primer artefact
                n_palin[jid] += 1
                continue
            cands[jid].append((f"{mode}:{i}", l, r, res[f"PRIMER_PAIR_{i}_PRODUCT_SIZE"],
                               round(res[f"PRIMER_PAIR_{i}_PENALTY"], 3), round(res[f"PRIMER_LEFT_{i}_TM"], 1),
                               round(res[f"PRIMER_RIGHT_{i}_TM"], 1)))
        if mode == "flank":
            res_flank = res
    res = res_flank
    explain[jid] = "; ".join(k.split("_")[1].lower() + ": " + str(res.get(k, "")) for k in
                              ("PRIMER_LEFT_EXPLAIN", "PRIMER_RIGHT_EXPLAIN", "PRIMER_PAIR_EXPLAIN"))
    print(jid, "candidate pairs:", len(cands[jid]), explain[jid], flush=True)

# exact genome-wide counts (both strands): grep -F over a one-line-per-chromosome copy of the reference, one
# query per grep (see count_one). (Overlapping self-occurrences of one pattern are counted once by grep -o; Primer3 excludes
# poly-X > 4, so this cannot turn a multi-hit primer into a single hit except in low-complexity sequence.)
import subprocess
import collections
ONELINE = os.environ.get("P44_GENOME_ONELINE",
                         "/tmp/claude-1001/-home-gao--claude/5253fed0-e021-4d8b-b988-67860909f347/scratchpad/GRCh38_oneline.txt")
if not os.path.exists(ONELINE):
    print("building one-line genome copy:", ONELINE, flush=True)
    with open(L.REF) as fi, open(ONELINE + ".part", "w") as fo:
        first = True
        for line in fi:
            if line.startswith(">"):
                fo.write(("" if first else "\n") + line[1:].split()[0] + "\t")
                first = False
            else:
                fo.write(line.rstrip().upper())
        fo.write("\n")
    os.rename(ONELINE + ".part", ONELINE)
full = sorted({p for v in cands.values() for _, l, r, *_ in v for p in (l, r)})
end15 = sorted({p[-15:] for p in full})
# ONE grep per query (forward + reverse complement): a multi-pattern grep -o lets an overlapping pattern consume a
# match and silently under-counts (found 2026-10-25: a primer shifted by 1 base masked its neighbour -> 0 hits).
from concurrent.futures import ThreadPoolExecutor


def count_one(q):
    pats = f"-e {q}" + (f" -e {rc(q)}" if rc(q) != q else "")
    out = subprocess.run(f"grep -o -F {pats} {ONELINE} | wc -l", shell=True, capture_output=True, text=True,
                         executable="/bin/bash").stdout.strip()
    return q, int(out or 0)


counts = collections.Counter()
with ThreadPoolExecutor(max_workers=8) as ex:
    for q, n in ex.map(count_one, sorted(set(full) | set(end15))):
        counts[q] = n
print("counted", len(counts), "queries", flush=True)

rows_c, rows_s = [], []
for e in events:
    jid = e["junction_id"]
    chosen = None
    for i, l, r, size, pen, tl, tr in cands[jid]:
        strict = all(counts[p] <= 1 and counts[p[-15:]] <= 1 for p in (l, r))
        # tier 2: one primer cannot prime the reference at all (crosses the junction: full 0, 3' 15-mer <= 1) and the
        # other has exactly one full-length reference site; off-target 3'-end matches of that partner alone give only
        # linear amplification, so an exponential product requires the rearranged allele
        jspec = [p for p in (l, r) if counts[p] == 0 and counts[p[-15:]] <= 1]
        tier2 = (not strict) and len(jspec) >= 1 and all(counts[p] == 1 for p in (l, r) if p not in jspec)
        status = "pass_strict" if strict else ("pass_junction_specific" if tier2 else "fail")
        rows_c.append((jid, i, l, r, size, pen, tl, tr, counts[l], counts[l[-15:]], counts[r], counts[r[-15:]], status))
        if status != "fail":
            rank = 0 if strict else 1
            if chosen is None or rank < chosen[-1]:
                chosen = (i, l, r, size, pen, tl, tr, status, rank)
    if chosen:
        i, l, r, size, pen, tl, tr, status, _ = chosen
        rows_s.append((jid, e["type"], e["breakend_A"], e["breakend_B"], e["samples_present"], "pair_selected_" + status.replace("pass_", ""),
                       l, tl, r, tr, size, pen,
                       f"design pass {i.split(':')[0]}; reference hits full/3'15-mer: left {counts[l]}/{counts[l[-15:]]}, "
                       f"right {counts[r]}/{counts[r[-15:]]}" + ("; one primer crosses the junction and cannot prime the "
                       "reference" if status == "pass_junction_specific" else "") +
                       ("; caution: one primer is mostly (TTAGGG)n/(CCCTAA)n and will anneal at every telomere; only its "
                        "3' end is chr11 sequence, so extension requires the junction, but background is expected"
                        if any(("TTAGGG" * 2) in p or ("CCCTAA" * 2) in p for p in (l, r)) else "")))
    else:
        why = []
        if "telomeric" in e.get("seqclass_flank_B", "") or "telomeric" in e.get("seqclass_flank_A", ""):
            why.append("one junction side is (TTAGGG)n telomeric repeat, present at every chromosome end; a "
                       "chr11-side primer combined with a telomere-repeat primer (single-telomere-length-type assay) "
                       "is the appropriate design")
        if not cands[jid] and n_palin[jid]:
            why.append(f"the junction template is an inverted duplication (fold-back): all {n_palin[jid]} Primer3 pairs had "
                       "primers that are reverse complements of each other, so standard two-primer junction PCR is not "
                       "applicable; an assay across the fold-back needs a different design (e.g. one primer outside the "
                       "duplicated segment combined with restriction digestion before PCR)")
        elif not cands[jid]:
            why.append("Primer3 returned no pair within the constraints (" + explain[jid] + ")")
        elif not why:
            why.append("in every candidate pair at least one primer has a full-length or 3' 15-mer match at more than one reference site (repetitive flank)")
        rows_s.append((jid, e["type"], e["breakend_A"], e["breakend_B"], e["samples_present"], "no_unique_pair",
                       "", "", "", "", "", "", "; ".join(why)))
with open(f"{OUT}/primer_candidates.tsv", "w") as fh:
    fh.write("event_id\tdesign_pass:rank\tleft_primer\tright_primer\tproduct_bp\tpair_penalty\tleft_tm\tright_tm\t"
             "left_full_hits\tleft_3prime15_hits\tright_full_hits\tright_3prime15_hits\tspecificity\n")
    for r in rows_c:
        fh.write("\t".join(map(str, r)) + "\n")
with open(f"{OUT}/primer_design.tsv", "w") as fh:
    fh.write("event_id\ttype\tbreakend_A\tbreakend_B\tsamples_present\tdesign_status\tleft_primer\tleft_tm\t"
             "right_primer\tright_tm\tproduct_bp_on_junction_template\tpair_penalty\tnote\n")
    for r in rows_s:
        fh.write("\t".join(map(str, r)) + "\n")
print("stage14 OK", [(r[0], r[5]) for r in rows_s])
