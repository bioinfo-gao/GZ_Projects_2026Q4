"""Stage G2b — genome-wide annotated per-event table for the DLD1 clone comparison (client Analysis 2 / deliverable 3).

Input: G2a dld1_S01_S05.somatic_tagged.bcf (every merged DELLY site re-genotyped in S01-S05; FILTER = DELLY somatic
filter tag, records never removed), P41 per-sample Manta calls (diploidSV any FILTER + candidateSV), GENCODE v45 genes.

Per event x sample (read counts are DELLY's: pairs DV/DR and junction reads RV/RR; not de-duplicated fragments):
  alt = DV + RV ; ref = DR + RR ; AF = alt / (alt + ref)
  present       genotype non-reference AND alt >= 3
  low_support   alt 1-2 (reported, never used to call an event present or absent)
  absent        alt == 0 AND ref >= 10   (AF 0.25 missed with probability <= 0.75^10 = 6%; same rule as the locus table)
  underpowered  alt == 0 AND ref < 10
Comparison (cases S04, S05; baseline S01, S02, S03):
  dox_specific               present in >= 1 case, absent in all 3 baseline samples
  dox_specific_underpowered  present in >= 1 case, no baseline with alt > 0, >= 1 baseline underpowered
  shared_with_baseline       present in >= 1 case and alt > 0 in >= 1 baseline sample
  baseline_only              present in >= 1 baseline sample only
  not_present                present in no sample
  present_in_both_dox        1 if present in S04 and S05
Manta column: samples whose own Manta calls match the event at both ends (+-100 bp, same orientation; svmatch.classify).
Annotation: GENCODE v45 genes (gene_type protein_coding or lncRNA) whose body contains each breakend.
Output: analysis/20261025_p45/G2_dld1_genome_wide/{dld1_genome_wide_events.tsv, dld1_dox_specific_events.tsv,
        summary_counts.tsv}
"""
import bisect
import collections
import gzip
import os
import sys
import pysam
sys.path.insert(0, os.path.dirname(__file__))
import p45lib as L
import svmatch as SM
from svmatch import pairs, classify

SM.TOL = 100
D = f"{L.OUT}/G2_dld1_genome_wide"
P41 = "/home/gao/projects_2026M7/41_Vedict_18_wgs/analysis/sv"
GTF = "/Work_bio/references/Homo_sapiens/GRCh38/human_gencode_v45/gencode.v45.annotation.gtf"
CASE = ["S04_DLD1_c3_dox", "S05_DLD1_c4_dox"]
BASE = ["S01_DLD1_pop", "S02_DLD1_c1_ctrl", "S03_DLD1_c2_ctrl"]
DLD1 = BASE + CASE
SH = {s: s.split("_")[0] for s in DLD1}
ALT_MIN, REF_MIN = 3, 10
pysam.set_verbosity(0)

# ---- genes
genes = collections.defaultdict(list)
with open(GTF) as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        f = line.split("\t", 9)
        if f[2] != "gene":
            continue
        gt = f[8].split('gene_type "')[1].split('"')[0]
        if gt not in ("protein_coding", "lncRNA"):
            continue
        genes[f[0]].append((int(f[3]), int(f[4]), f[8].split('gene_name "')[1].split('"')[0]))
gstart = {c: [g[0] for g in sorted(v)] for c, v in genes.items()}
for c in genes:
    genes[c].sort()


def genes_at(c, p):
    out = []
    i = bisect.bisect_right(gstart.get(c, []), p)
    for s, e, n in genes.get(c, [])[max(0, i - 400):i]:      # genes start <= p; longest human gene < 2.5 Mb
        if s <= p <= e:
            out.append(n)
    return ",".join(sorted(set(out)))


# ---- Manta per-sample
manta = {s: [pysam.VariantFile(f"{P41}/manta_gw/{s}/results/variants/{f}") for f in ("diploidSV.vcf.gz", "candidateSV.vcf.gz")]
         for s in DLD1}


def manta_match(s, A, B):
    recs = []
    for vf in manta[s]:
        for b in (A, B):
            try:
                for rec in vf.fetch(b[0], max(0, b[1] - SM.TOL - 1), b[1] + SM.TOL):
                    recs += pairs(rec)
            except (ValueError, KeyError):
                pass
    return classify((A, B), recs)[0] == "two_end"


def state(gt, alt, ref):
    nonref = gt is not None and any(g not in (0, None) for g in gt)
    if alt >= ALT_MIN and nonref:
        return "present"
    if alt >= 1:
        return "low_support"
    return "absent" if ref >= REF_MIN else "underpowered"


vf = pysam.VariantFile(f"{D}/dld1_S01_S05.somatic_tagged.bcf")
assert set(vf.header.samples) == set(DLD1), list(vf.header.samples)
rows = []
for rec in vf:
    pr = pairs(rec)
    if not pr or pr[0][0] in ("INS", "SGL"):
        A = B = (rec.chrom, rec.pos, "+")
        ins = True
    else:
        A, B = pr[0]
        ins = False
    st, af, cnt = {}, {}, {}
    for s in DLD1:
        f = rec.samples[s]
        alt = (f.get("DV") or 0) + (f.get("RV") or 0)
        ref = (f.get("DR") or 0) + (f.get("RR") or 0)
        st[s] = state(f.get("GT"), alt, ref)
        cnt[s] = (alt, ref)
        af[s] = f"{alt / (alt + ref):.3f}" if alt + ref else "NA"
    cp = [s for s in CASE if st[s] == "present"]
    bp = [s for s in BASE if st[s] == "present"]
    b_any_alt = [s for s in BASE if cnt[s][0] > 0]
    if cp and not b_any_alt:
        lab = "dox_specific" if all(st[s] == "absent" for s in BASE) else "dox_specific_underpowered"
    elif cp:
        lab = "shared_with_baseline"
    elif bp:
        lab = "baseline_only"
    else:
        lab = "not_present"
    msamp = [SH[s] for s in DLD1 if not ins and manta_match(s, A, B)]
    dsamp = [SH[s] for s in DLD1 if st[s] == "present"]
    callers = "+".join(x for x, ok in (("DELLY", bool(dsamp)), ("Manta", bool(msamp))) if ok) or "DELLY_low_support_only"
    info = rec.info
    row = collections.OrderedDict(
        event_id=rec.id, svtype=info.get("SVTYPE"), chrom_A=A[0], pos_A=A[1], orient_A=A[2], chrom_B=B[0], pos_B=B[1],
        orient_B=B[2], size_bp=(abs(B[1] - A[1]) if A[0] == B[0] and not ins else ""),
        precise=int("PRECISE" in info), site_quality=rec.qual, delly_somatic_filter=";".join(rec.filter.keys()),
        caller_support=callers, delly_present_samples=",".join(dsamp), manta_two_end_samples=",".join(msamp),
        comparison_S04_S05_vs_S01_S03=lab, present_in_both_dox=int(len(cp) == 2),
        genes_at_A=genes_at(A[0], A[1]), genes_at_B=genes_at(B[0], B[1]))
    for s in DLD1:
        row[f"{SH[s]}_alt"] = cnt[s][0]
        row[f"{SH[s]}_ref"] = cnt[s][1]
        row[f"{SH[s]}_AF"] = af[s]
        row[f"{SH[s]}_state"] = st[s]
    rows.append(row)

cols = list(rows[0].keys())
order = {"dox_specific": 0, "dox_specific_underpowered": 1, "shared_with_baseline": 2, "baseline_only": 3, "not_present": 4}
rows.sort(key=lambda r: (order[r["comparison_S04_S05_vs_S01_S03"]], L.chrom_rank(r["chrom_A"]), r["pos_A"]))
with open(f"{D}/dld1_genome_wide_events.tsv", "w") as fh:
    fh.write("\t".join(cols) + "\n")
    for r in rows:
        fh.write("\t".join(map(str, r.values())) + "\n")
sel = [r for r in rows if r["comparison_S04_S05_vs_S01_S03"].startswith("dox_specific")]
with open(f"{D}/dld1_dox_specific_events.tsv", "w") as fh:
    fh.write("\t".join(cols) + "\n")
    for r in sel:
        fh.write("\t".join(map(str, r.values())) + "\n")
summ = [("events_total", len(rows))]
summ += [(f"comparison_{k}", v) for k, v in collections.Counter(r["comparison_S04_S05_vs_S01_S03"] for r in rows).items()]
summ += [("dox_specific_both_dox", sum(r["present_in_both_dox"] for r in sel if r["comparison_S04_S05_vs_S01_S03"] == "dox_specific"))]
summ += [(f"dox_specific_svtype_{k}", v) for k, v in collections.Counter(r["svtype"] for r in sel if r["comparison_S04_S05_vs_S01_S03"] == "dox_specific").items()]
summ += [("dox_specific_PASS_somatic_filter", sum(1 for r in sel if r["comparison_S04_S05_vs_S01_S03"] == "dox_specific" and r["delly_somatic_filter"] == "PASS"))]
summ += [("dox_specific_with_Manta", sum(1 for r in sel if r["comparison_S04_S05_vs_S01_S03"] == "dox_specific" and r["manta_two_end_samples"]))]
summ += [("delly_somatic_PASS_total", sum(1 for r in rows if r["delly_somatic_filter"] == "PASS"))]
with open(f"{D}/summary_counts.tsv", "w") as fh:
    fh.write("metric\tvalue\n")
    for k, v in summ:
        fh.write(f"{k}\t{v}\n")
print("\n".join(f"{k}\t{v}" for k, v in summ))
print("G2b OK")
