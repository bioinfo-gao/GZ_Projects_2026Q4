"""Stage 9b-2 — chr11 B-allele frequency (BAF) and loss of heterozygosity (LOH) for S01-S10.

Informative sites: biallelic SNVs (QUAL >= 30) that are heterozygous in EVERY baseline sample of the cell line
(DLD1: S01, S02, S03; GM12878: S06), each baseline with depth >= 10 and alt fraction 0.25-0.75.
For each sample at those sites: BAF = AD_alt / (AD_ref + AD_alt) (depth >= 8).
Segments (as P41 section 7.2): chr11p 0-50 Mb (internal control), proximal 99.0-114.5 Mb, distal 114.7-135.0 Mb.
Per segment: n sites, median |BAF - 0.5| (allelic imbalance), fraction of sites with BAF < 0.1 or > 0.9 (LOH-like).
Expectation if one of two copies is lost distal to the step: distal sites become LOH-like in S04/S05.
Outputs: analysis/<RUN>/stage9_supporting/chr11_baf/{baf_sites.tsv.gz, baf_segments.tsv}
"""
import gzip
import os
import sys
import numpy as np
import pysam
sys.path.insert(0, os.path.dirname(__file__))
import p44lib as L

RUN = os.environ.get("P44_RUN_DATE", "20261024")
BASE = f"{L.PROJ}/analysis/{RUN}_read_rescue/stage9_supporting"
OUT = f"{BASE}/chr11_baf"
os.makedirs(OUT, exist_ok=True)
VCF = f"{BASE}/chr11_snv/chr11_S01_S10.QUAL30.vcf.gz"
SEGS = {"chr11p_0-50Mb": (0, 50_000_000), "proximal_99.0-114.5Mb": (99_000_000, 114_500_000),
        "distal_114.7-135.0Mb": (114_700_000, 135_000_000)}
LINES = {"DLD1-136": (["S01_DLD1_pop", "S02_DLD1_c1_ctrl", "S03_DLD1_c2_ctrl"],
                      ["S01_DLD1_pop", "S02_DLD1_c1_ctrl", "S03_DLD1_c2_ctrl", "S04_DLD1_c3_dox", "S05_DLD1_c4_dox"]),
         "GM12878": (["S06_GM12878_pop"],
                     ["S06_GM12878_pop", "S07_GM12878_c5", "S08_GM12878_c6", "S09_GM12878_c7", "S10_GM12878_c8"])}

vf = pysam.VariantFile(VCF)
data = {line: [] for line in LINES}
for rec in vf:
    if len(rec.alts or ()) != 1 or len(rec.ref) != 1 or len(rec.alts[0]) != 1:
        continue
    for line, (base, tests) in LINES.items():
        ok = True
        for s in base:
            ad = rec.samples[s].get("AD")
            if ad is None or None in ad or sum(ad) < 10 or not (0.25 <= ad[1] / sum(ad) <= 0.75):
                ok = False
                break
        if not ok:
            continue
        row = [rec.pos]
        for s in tests:
            ad = rec.samples[s].get("AD")
            row.append(ad[1] / sum(ad) if ad and None not in ad and sum(ad) >= 8 else np.nan)
        data[line].append(row)

with gzip.open(f"{OUT}/baf_sites.tsv.gz", "wt") as fh, open(f"{OUT}/baf_segments.tsv", "w") as fs:
    fs.write("cell_line\tsample\tsegment\tn_sites\tmedian_abs_BAF_minus_0.5\tfrac_LOH_like\n")
    for line, (base, tests) in LINES.items():
        a = np.array(data[line], dtype=float)
        fh.write(f"#{line}\tpos\t" + "\t".join(tests) + "\n")
        for r in a:
            fh.write(f"{line}\t{int(r[0])}\t" + "\t".join("" if np.isnan(x) else f"{x:.3f}" for x in r[1:]) + "\n")
        for seg, (s0, s1) in SEGS.items():
            m = (a[:, 0] >= s0) & (a[:, 0] < s1)
            for k, s in enumerate(tests, 1):
                v = a[m, k]
                v = v[~np.isnan(v)]
                fs.write(f"{line}\t{s}\t{seg}\t{len(v)}\t{np.median(np.abs(v - 0.5)) if len(v) else float('nan'):.3f}\t"
                         f"{np.mean((v < 0.1) | (v > 0.9)) if len(v) else float('nan'):.3f}\n")
print(open(f"{OUT}/baf_segments.tsv").read())
print("stage9b-2 OK")

# ---- which homolog was retained distal to the step: S04 vs S05 retained-allele concordance (added 2026-10-25)
a = np.array(data["DLD1-136"], dtype=float)
m = (a[:, 0] >= 114_700_000) & (a[:, 0] < 135_000_000)
b4, b5 = a[m, 4], a[m, 5]                        # columns: pos, S01, S02, S03, S04, S05
ok = ~np.isnan(b4) & ~np.isnan(b5) & ((b4 < 0.1) | (b4 > 0.9)) & ((b5 < 0.1) | (b5 > 0.9))
same = ((b4[ok] > 0.9) == (b5[ok] > 0.9)).sum()
with open(f"{OUT}/retained_homolog_S04_vs_S05.tsv", "w") as fh:
    fh.write("segment\tsites_LOH_in_both\tsame_allele_retained\tfraction_same\n")
    fh.write(f"distal_114.7-135.0Mb\t{ok.sum()}\t{same}\t{same / max(ok.sum(), 1):.4f}\n")
print(open(f"{OUT}/retained_homolog_S04_vs_S05.tsv").read())
