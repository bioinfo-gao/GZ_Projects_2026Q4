"""Stage 10 — final tables (client deliverables) incl. supporting Analysis 2 part 9a (clone-specific labels).

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
import p44lib as L

RUN = os.environ.get("P44_RUN_DATE", "20261024")
BASE = f"{L.PROJ}/analysis/{RUN}_read_rescue"
OUT = f"{BASE}/stage10_tables{L.SUFFIX}"
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
conc = {d["junction_id"]: d for d in tsv(f"{BASE}/stage8_concordance{L.SUFFIX}/concordance_summary.tsv")}
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
        junction_id=jid, type=a["type"], breakend_A=a["breakend_A"], breakend_B=a["breakend_B"],
        size_bp=a["size_bp"], classification=a["classification"], ambiguity_reason=a["ambiguity_reason"],
        n_alternative_representations=a["n_alternative_representations"],
        alternative_representations=a["alternative_representations"],
        seqclass_flank_A=a["seqclass_flank_A"], seqclass_flank_B=a["seqclass_flank_B"],
        caller_support=c.get("caller_support", "NA"),
        manta_samples=c.get("manta_unfiltered_or_joint_samples", ""), manta_PASS_samples=c.get("manta_pass_samples", ""),
        delly_samples=c.get("delly_unfiltered_or_regeno_samples", ""), delly_PASS_samples=c.get("delly_pass_samples", ""),
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
        P41_event=a["P41_event"], P41_max_tier=a["P41_max_tier"])
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
with open(f"{OUT}/summary_counts.tsv", "w") as fh:
    fh.write("metric\tvalue\n")
    for k, v in summ:
        fh.write(f"{k}\t{v}\n")
print("\n".join(f"{k}\t{v}" for k, v in summ))
print("stage10 OK")
