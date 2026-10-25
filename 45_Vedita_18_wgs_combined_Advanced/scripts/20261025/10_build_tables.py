"""Stage 10 (P45) — final tables. Copied from P44 10_build_tables.py (from_P44/); P45 changes, marked "P45:" below:
  * inputs are the inherited P44 round-2 products (analysis/20261025_inherited/P44_r2/);
  * caller columns come from Stage 8x (Manta, DELLY, GRIDSS, SvABA) and claim_level from Stage 9x;
  * client-facing junction IDs drop the internal prefix (P44J0332 -> J0332); columns named after internal project
    numbers are renamed after the 2026-10-09 report the client holds;
  * the read ledger adds strand, read-orientation query interval, SEQ, a per-event fragment ID and, per read end,
    the query overlap (+, microhomology) or gap (-, inserted bases) between consecutive split segments.

Original P44 description:

Inputs: Stage 7 junction_annotated.tsv, Stage 8 concordance_summary.tsv, Stage 6 junction_by_sample.tsv,
Stage 3 filter_audit.tsv, P42 fixed-candidate coordinates (comparison annotation only; P44 is independent).

9a comparison labels (criteria stated so the client can re-filter):
  present in a sample   >= 1 non-duplicate after-rescue fragment (any MAPQ)
  confident in a sample >= 2 such fragments
  absent in a sample    0 fragments AND >= 10 reference-spanning fragments at the in-target breakend
                        (allele fraction 0.25 missed with probability <= 0.75^10 = 6%)
  underpowered          0 fragments and < 10 reference-spanning fragments
  DLD1 (+Dox S04,S05 vs S01,S02,S03) and GM12878 (each clone S07-S10 vs S06):
    case_specific            present in >= 1 case sample, absent in ALL baseline samples
    case_specific_underpowered present in >= 1 case sample, no baseline present, >= 1 baseline underpowered
    shared_with_baseline     present in >= 1 case and >= 1 baseline sample
    baseline_only            present only in baseline sample(s)
    not_present              present in no sample of the cell line
  No allele-fraction minimum is applied; allele fraction and its 95% CI are reported per sample.
Outputs: analysis/<RUN>/stage10_tables{L.SUFFIX}/
  junction_table_target.tsv        one row per junction with an end in the target region (main table)
  junction_table_window_context.tsv junctions with ends only in the W flanks (outside target)
  junction_by_sample_long.tsv      per junction x sample counts (S01-S10 + background S11-S18)
  rescue_source_summary.tsv        per sample: QNAMEs by rescue source; completeness
  summary_counts.tsv               headline numbers used in the report (single source for report text)
"""
import collections
import gzip
import os
import sys
sys.path.insert(0, os.path.dirname(__file__))
import p45lib as L

BASE = L.P44R2                                   # P45: inherited P44 round-2 stage 1-8
OUT = f"{L.OUT}/stage10_tables"
os.makedirs(OUT, exist_ok=True)
ALL = [r["sample_id"] for r in L.samples()]
REP = [r["sample_id"] for r in L.samples(reported_only=True)]
SHORT = {s: s.split("_")[0] for s in ALL}


def tsv(path):
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt") as fh:
        h = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            yield dict(zip(h, line.rstrip("\n").split("\t")))


ann = {d["junction_id"]: d for d in tsv(f"{BASE}/stage7_annotation{L.SUFFIX}/junction_annotated.tsv")}
conc = {d["junction_id"]: d for d in tsv(f"{L.OUT}/stage8x_concordance/concordance_summary_4callers.tsv")}   # P45
claim = {d["junction_id"]: d for d in tsv(f"{L.OUT}/stage9x_coherence/claim_ladder.tsv")}                    # P45
by = collections.defaultdict(dict)
for d in tsv(f"{BASE}/stage6_genotype{L.SUFFIX}/junction_by_sample.tsv"):
    by[d["junction_id"]][d["sample"]] = d


def state(d):
    f = int(d["after_frag"])
    if f >= 1:
        return "present"
    return "absent" if d["status"] == "absent" else "underpowered"


def compare(jid, cases, bases):
    cases = [s for s in cases if s in by[jid]]          # test mode may run on a sample subset
    bases = [s for s in bases if s in by[jid]]
    if not cases or not bases:
        return "not_evaluated", ""
    st = {s: state(by[jid][s]) for s in cases + bases}
    cp = [s for s in cases if st[s] == "present"]
    bp = [s for s in bases if st[s] == "present"]
    if cp and not bp:
        lab = "case_specific" if all(st[s] == "absent" for s in bases) else "case_specific_underpowered"
    elif cp and bp:
        lab = "shared_with_baseline"
    elif bp:
        lab = "baseline_only"
    else:
        lab = "not_present"
    return lab, ",".join(SHORT[s] for s in cp)


rows_t, rows_w = [], []
for jid, a in ann.items():
    c = conc.get(jid, {})
    lab_d, case_d = compare(jid, ["S04_DLD1_c3_dox", "S05_DLD1_c4_dox"],
                            ["S01_DLD1_pop", "S02_DLD1_c1_ctrl", "S03_DLD1_c2_ctrl"])
    gm = []
    for s in ["S07_GM12878_c5", "S08_GM12878_c6", "S09_GM12878_c7", "S10_GM12878_c8"]:
        lab, _ = compare(jid, [s], ["S06_GM12878_pop"])
        if lab.startswith("case_specific"):
            gm.append(f"{SHORT[s]}:{lab}")
    both_dox = all(int(by[jid][s]["after_frag"]) >= 1 for s in ["S04_DLD1_c3_dox", "S05_DLD1_c4_dox"] if s in by[jid])
    row = collections.OrderedDict(
        junction_id=L.cid(jid), type=a["type"], breakend_A=a["breakend_A"], breakend_B=a["breakend_B"],
        size_bp=a["size_bp"], classification=a["classification"], claim_level=claim[jid]["claim_level"],
        structurally_coherent_in=claim[jid]["coherent_in_samples"], ambiguity_reason=a["ambiguity_reason"],
        n_alternative_representations=a["n_alternative_representations"],
        alternative_representations=L.cid(a["alternative_representations"]),
        seqclass_flank_A=a["seqclass_flank_A"], seqclass_flank_B=a["seqclass_flank_B"],
        caller_support=c.get("caller_support_4", "NA"), n_callers=c.get("n_callers_two_end", "NA"),
        manta_samples=c.get("manta_samples", ""), manta_PASS_samples=c.get("manta_PASS_samples", ""),
        delly_samples=c.get("delly_samples", ""), delly_PASS_samples=c.get("delly_PASS_samples", ""),
        gridss_samples=c.get("gridss_samples", ""), gridss_PASS_samples=c.get("gridss_PASS_samples", ""),
        svaba_samples_S01_S05=c.get("svaba_samples_S01_S05", ""),
        svaba_PASS_samples_S01_S05=c.get("svaba_PASS_samples_S01_S05", ""),
        caller_one_end_anchor_samples=c.get("any_one_end_anchor", ""),
        evidence_types=a["evidence_types"], max_after_frag_S01_S10=a["max_after_frag_S01_S10"],
        samples_present=a["samples_present_S01_S10"], DLD1_dox_comparison=lab_d, DLD1_dox_case_samples=case_d,
        present_in_both_dox_clones=int(both_dox), GM12878_clone_specific=";".join(gm),
        background_S11_S18_samples_present=a["background_samples_present_S11_S18"],
        cell_lines_present=a["cell_lines_present"],
        microhomology_bp=a["microhomology_bp"], microhomology_seq=a["microhomology_seq"],
        inserted_bp=a["inserted_bp"], inserted_seq=a["inserted_seq"],
        breakpoint_interval_A=a["uncertainty_A"], breakpoint_interval_B=a["uncertainty_B"],
        flankA_n_equiv_placements=a["flankA_n_equiv"], flankB_n_equiv_placements=a["flankB_n_equiv"],
        flankA_competing_placements=a["flankA_competitors"], flankB_competing_placements=a["flankB_competitors"],
        repeat_A=a["repeat_A"], repeat_B=a["repeat_B"], segdup_A=a["segdup_A"], segdup_B=a["segdup_B"],
        simple_repeat_A=a["simple_A"], simple_repeat_B=a["simple_B"],
        breakpoint_source=a["rep_source"], n_assembled_contigs=a["n_contigs"],
        event_in_20261009_report=a["P41_event"], tier_in_20261009_report=a["P41_max_tier"])
    for s in REP:
        if s not in by[jid]:
            continue
        d = by[jid][s]
        row[f"{SHORT[s]}_before"] = d["before_frag"]
        row[f"{SHORT[s]}_after"] = d["after_frag"]
        row[f"{SHORT[s]}_ref"] = d["refA_frag"] if L.in_w(*[a["breakend_A"].rsplit(":", 2)[0],
                                                           int(a["breakend_A"].rsplit(":", 2)[1])]) else d["refB_frag"]
        row[f"{SHORT[s]}_AF"] = d["af_W_side"]
        row[f"{SHORT[s]}_AF95CI"] = f"{d['af_lo95']}-{d['af_hi95']}"
        row[f"{SHORT[s]}_state"] = state(d)
    (rows_t if a["in_target"] == "1" else rows_w).append(row)

order = {"Well_supported": 0, "Limited": 1, "Ambiguous_repeat": 2, "Weak": 3, "Not_supported": 4}
for rows, fn in ((rows_t, "junction_table_target.tsv"), (rows_w, "junction_table_window_context.tsv")):
    rows.sort(key=lambda r: (order[r["classification"]], -int(r["max_after_frag_S01_S10"]), r["junction_id"]))
    with open(f"{OUT}/{fn}", "w") as fh:
        if rows:
            fh.write("\t".join(rows[0].keys()) + "\n")
            for r in rows:
                fh.write("\t".join(map(str, r.values())) + "\n")

with open(f"{OUT}/junction_by_sample_long.tsv", "w") as fh:
    cols = None
    for jid in ann:
        for s in ALL:
            d = dict(by[jid][s])
            d["junction_id"] = L.cid(d["junction_id"])                       # P45
            d["event_group"] = L.cid(d.get("event_group", ""))
            d["panel"] = "reported" if s in REP else "background"
            if cols is None:
                cols = list(d.keys())
                fh.write("\t".join(cols) + "\n")
            fh.write("\t".join(str(d[k]) for k in cols) + "\n")

src_rows = []
for s in ALL:
    a = {(d["step"], d["item"]): d["count"] for d in tsv(f"{BASE}/stage3_rescue/{s}/filter_audit.tsv")}
    src_rows.append([s, a[("input_records", "window")], a[("qnames", "total")]] +
                    [a[("qnames_by_source", x)] for x in "abcde"] + [a[("qnames_only_source", x)] for x in "abcde"] +
                    [a[("qnames", "complete_pairs_and_SA")], a[("qnames", "incomplete")],
                     a[("unmapped_pairs", "reads_screened")]])
with open(f"{OUT}/rescue_source_summary.tsv", "w") as fh:
    fh.write("sample\trecords_in_W\tqnames_total\t" + "\t".join(f"qnames_with_{x}" for x in "abcde") + "\t" +
             "\t".join(f"qnames_only_{x}" for x in "abcde") + "\tqnames_complete\tqnames_incomplete\tunmapped_pair_reads_screened\n")
    for r in src_rows:
        fh.write("\t".join(map(str, r)) + "\n")

# headline numbers (single source for report text)
cnt = collections.Counter(r["classification"] for r in rows_t)
summ = [("junctions_in_target_total", len(rows_t)), ("junctions_window_context_total", len(rows_w))]
summ += [(f"target_{k}", cnt.get(k, 0)) for k in order]
for s in REP:
    summ.append((f"{SHORT[s]}_target_junctions_before_ge1", sum(int(r[f"{SHORT[s]}_before"]) >= 1 for r in rows_t)))
    summ.append((f"{SHORT[s]}_target_junctions_after_ge1", sum(int(r[f"{SHORT[s]}_after"]) >= 1 for r in rows_t)))
    summ.append((f"{SHORT[s]}_target_junctions_after_ge2", sum(int(r[f"{SHORT[s]}_after"]) >= 2 for r in rows_t)))
summ += [(f"target_caller_{k}", v) for k, v in collections.Counter(r["caller_support"] for r in rows_t).items()]
summ += [(f"target_DLD1_{k}", v) for k, v in collections.Counter(r["DLD1_dox_comparison"] for r in rows_t).items()]
summ.append(("target_present_in_both_dox", sum(int(r["present_in_both_dox_clones"]) for r in rows_t)))
summ += [(f"target_claim_{k}", v) for k, v in collections.Counter(r["claim_level"] for r in rows_t).items()]   # P45
summ += [(f"window_claim_{k}", v) for k, v in collections.Counter(r["claim_level"] for r in rows_w).items()]
with open(f"{OUT}/summary_counts.tsv", "w") as fh:
    fh.write("metric\tvalue\n")
    for k, v in summ:
        fh.write(f"{k}\t{v}\n")
print("\n".join(f"{k}\t{v}" for k, v in summ))
print("stage10 OK")

# ------------------------------------------------------------------ P3: read-level evidence ledger for key events
# Every event with >= 2 fragments in a reported sample: every supporting read end / pair (Stage 6 support, incl.
# alternative partner placements of the event) with EVERY original alignment record of that read name.
import pysam
SUFFIX_STAGE6 = f"{BASE}/stage6_genotype{L.SUFFIX}"
members_of = collections.defaultdict(list)
for d in tsv(f"{SUFFIX_STAGE6}/junctions.tsv"):
    members_of[d["event_group"]].append(d["junction_id"])
key_events = [dict(r, junction_id=next(j for j in ann if L.cid(j) == r["junction_id"]))      # P45: back to internal id
              for r in rows_t + rows_w if int(r["max_after_frag_S01_S10"]) >= 2]
want = collections.defaultdict(list)          # (sample, qname) -> [(event, read_end, class, AS_t, AS_g)]
for d in tsv(f"{SUFFIX_STAGE6}/read_support.tsv.gz"):
    if d["class"] not in ("alt", "alt_disc"):
        continue
    for ev in key_events:
        if d["junction_id"] in members_of[ev["junction_id"]]:
            want[(d["sample"], d["qname"])].append((ev["junction_id"], d["read_end"], d["class"], d["AS_template"],
                                                   d["AS_genome"], d["junction_id"]))
# one entry per (event, read end, class): the representation on which that read end scored best
for k, v in want.items():
    best = {}
    for e in v:
        kk = (e[0], e[1], e[2])
        sc = int(e[3]) if e[3] not in ("", None) else -1
        if kk not in best or sc > best[kk][0]:
            best[kk] = (sc, e)
    want[k] = [b[1] for b in best.values()]
ledger = []
for s in sorted({k[0] for k in want}):
    qs = {k[1] for k in want if k[0] == s}
    for r in pysam.AlignmentFile(f"{BASE}/stage3_rescue/{s}/rescue.bam").fetch(until_eof=True):
        if r.query_name not in qs:
            continue
        for ev, re_, cls, ast, asg, via in want[(s, r.query_name)]:
            # P45: strand, query interval in read orientation, SEQ (primary/supplementary carry it; secondary '*')
            if r.is_unmapped:
                q0 = q1 = ""
            else:
                g = L.seg_from_read(r)
                q0, q1 = g.qs, g.qe
            ledger.append((L.cid(ev), L.cid(via), s, r.query_name, cls, re_, 1 if r.is_read1 else 2, r.flag,
                           "secondary" if r.is_secondary else "supplementary" if r.is_supplementary else
                           ("unmapped" if r.is_unmapped else "primary"),
                           r.reference_name if not r.is_unmapped else "*", r.reference_start + 1, r.cigarstring or "*",
                           r.mapping_quality, int(r.is_duplicate), r.get_tag("SA") if r.has_tag("SA") else "",
                           r.get_tag("XA") if r.has_tag("XA") else "", ast, asg,
                           r.get_tag("ZR") if r.has_tag("ZR") else "",
                           "*" if r.is_unmapped else ("-" if r.is_reverse else "+"), q0, q1,
                           r.query_sequence or "*"))
ledger.sort(key=lambda x: (x[0], x[2], x[3], x[6], x[10]))
# P45: fragment ID per event (F01, F02 ... in sample/read-name order) and, per read end, the overlap (+) / gap (-)
# between consecutive primary/supplementary segments in read orientation (P42 C2 read-end resolution, generalised)
frag_id, segs = {}, collections.defaultdict(list)
for x in ledger:
    k = (x[0], x[2], x[3])
    if k not in frag_id:
        frag_id[k] = f"F{sum(1 for kk in frag_id if kk[0] == x[0]) + 1:02d}"
    if x[8] in ("primary", "supplementary") and x[20] != "":
        segs[(x[0], x[2], x[3], x[6])].append((x[20], x[21]))
ovl = {}
for k, v in segs.items():
    v = sorted(set(v))
    ovl[k] = ";".join(str(v[i][1] - v[i + 1][0]) for i in range(len(v) - 1)) if len(v) > 1 else ""
seen = set()
with open(f"{OUT}/key_junction_read_evidence.tsv", "w") as fh:
    fh.write("event_id\tsupporting_junction_representation\tsample\tfragment_id\tread_name\tsupport_class\t"
             "supporting_read_end\trecord_read_end\tflag\trecord_class\tchrom\tpos\tcigar\tmapq\tduplicate\tSA\tXA\t"
             "AS_on_junction_template\tAS_best_genome\trescue_sources\tstrand\tquery_start_read_orientation\t"
             "query_end_read_orientation\tsplit_query_overlap_bp_this_read_end\tsequence\n")
    for x in ledger:
        if x in seen:
            continue
        seen.add(x)
        y = list(x[:3]) + [frag_id[(x[0], x[2], x[3])]] + list(x[3:22]) + [ovl.get((x[0], x[2], x[3], x[6]), ""), x[22]]
        fh.write("\t".join(map(str, y)) + "\n")
print("key-event ledger:", len(seen), "rows,", len(key_events), "events")
