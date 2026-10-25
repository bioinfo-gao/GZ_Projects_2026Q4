"""Stage 10c (P45, INTERNAL ONLY — never copied into the client package) — P42 C1-C4 vs P44/P45 junctions.

For each P42 fixed candidate, the P45 junction whose two breakends both lie within 50 bp of the candidate's two
endpoints (any order; P42 endpoints carry no orientation) is reported next to P42's own numbers:
  P42: raw / standard-eligible dual-endpoint QNAMEs per DLD1 sample, Manta/DELLY exact counts, rescue impact
  P45: evidence class, claim level, before/after fragments per DLD1 sample, 4-caller support (Stage 8x, with distance)
If no P45 junction matches, the nearest P45 junction with one end within 1 kb is listed so the reason is visible.
Output: docs/P45_P42_P44_reconciliation.tsv  (+ printed summary for the operations log)
"""
import collections
import csv
import os
import sys
sys.path.insert(0, os.path.dirname(__file__))
import p45lib as L

R = f"{L.P42IN}/results"
DLD1 = ["S01", "S02", "S03", "S04", "S05"]


def tsv(p):
    with open(p) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def bk(s):
    c, p, _ = s.rsplit(":", 2)
    return c, int(p)


cands = {}
for d in tsv(f"{R}/candidate_manta_delly_concordance.tsv"):
    cands.setdefault(d["candidate_label"], [tuple((x.split(":")[0], int(x.split(":")[1]))) for x in d["fixed_endpoints"].split(";")])
p42c = collections.defaultdict(dict)
for d in tsv(f"{R}/clone_specific_comparison_with_criteria.tsv"):
    p42c[d["candidate_label"]][d["sample_id"][:3]] = (d["raw_dual_endpoint_fragment_count"], d["standard_eligible_dual_fragment_count"])
p42cc = collections.defaultdict(lambda: [0, 0])
for d in tsv(f"{R}/candidate_manta_delly_concordance.tsv"):
    p42cc[d["candidate_label"]][0] += int(d["manta_exact_two_end_record_count"])
    p42cc[d["candidate_label"]][1] += int(d["delly_exact_two_end_record_count"])
resc = {d["candidate_label"]: (d["raw_dual_endpoint_qname_count_before_rescue"], d["raw_dual_endpoint_qname_count_after_rescue"])
        for d in tsv(f"{R}/read_rescue_candidate_impact.tsv")}
ev = [r for f in ("junction_table_target.tsv", "junction_table_window_context.tsv") for r in tsv(f"{L.OUT}/stage10_tables/{f}")]
conc = collections.defaultdict(list)
for d in tsv(f"{L.OUT}/stage8x_concordance/concordance_long_4callers.tsv"):
    if d["match"] == "two_end" and d["sample"][:3] in DLD1:
        conc[L.cid(d["junction_id"])].append(f"{d['callset']}:{d['sample'][:3]}:{d['best_two_end_dist_bp']}bp")

rows = []
for lab, (e1, e2) in sorted(cands.items()):
    def match(r, tol):
        a, b = bk(r["breakend_A"]), bk(r["breakend_B"])
        ok = lambda x, y: x[0] == y[0] and abs(x[1] - y[1]) <= tol
        return (ok(a, e1) and ok(b, e2)) or (ok(a, e2) and ok(b, e1))
    hit = [r for r in ev if match(r, 50)]
    near = [] if hit else sorted((r for r in ev if any(bk(r[k])[0] == e[0] and abs(bk(r[k])[1] - e[1]) <= 1000
                                                    for k in ("breakend_A", "breakend_B") for e in (e1, e2))),
                                 key=lambda r: -int(r["max_after_frag_S01_S10"]))[:1]
    r = (hit or near or [None])[0]
    rows.append(collections.OrderedDict(
        P42_candidate=lab, P42_endpoints=f"{e1[0]}:{e1[1]};{e2[0]}:{e2[1]}",
        P42_raw_S01_S05=",".join(p42c[lab].get(s, ("NA", ""))[0] for s in DLD1),
        P42_standard_eligible_S01_S05=",".join(p42c[lab].get(s, ("", "NA"))[1] for s in DLD1),
        P42_manta_exact_total=p42cc[lab][0], P42_delly_exact_total=p42cc[lab][1],
        P42_rescue_before_after="/".join(resc.get(lab, ("NA", "NA"))),
        P45_match="two_end_within_50bp" if hit else ("nearest_one_end_within_1kb" if near else "none"),
        P45_junction=r["junction_id"] if r else "", P45_breakends=f"{r['breakend_A']};{r['breakend_B']}" if r else "",
        P45_class=r["classification"] if r else "", P45_claim=r["claim_level"] if r else "",
        P45_before_S01_S05=",".join(r[f"{s}_before"] for s in DLD1) if r else "",
        P45_after_S01_S05=",".join(r[f"{s}_after"] for s in DLD1) if r else "",
        P45_caller_support=r["caller_support"] if r else "",
        P45_two_end_caller_records_DLD1=";".join(conc.get(r["junction_id"], [])) if r else ""))
os.makedirs(f"{L.PROJ}/docs", exist_ok=True)
with open(f"{L.PROJ}/docs/P45_P42_P44_reconciliation.tsv", "w") as fh:
    fh.write("\t".join(rows[0].keys()) + "\n")
    for r in rows:
        fh.write("\t".join(map(str, r.values())) + "\n")
for r in rows:
    print(dict(r))
print("stage10c OK")
