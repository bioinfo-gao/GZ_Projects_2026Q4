"""Stage 9c — localise the S04/S05 depth step at the locus (supports Analysis 1: a junction that explains the
copy-number change must lie inside this interval).

Depth per base (samtools depth -a; default excludes unmapped/secondary/QC-fail/duplicate) at MAPQ >= 20 and
MAPQ >= 0 over chr11:114,580,000-114,650,000 for S01-S10, 200 bp bins, each divided by the sample's own
median over the flanking 114.40-114.58 Mb (local normalisation, same chromosome). Case ratio = case / mean of
the DLD1 no-dox samples (S01-S03); GM12878 clones vs S06 as negative comparison.
Step model: two-level fit on log2 ratio over 114.590-114.640 Mb; boundary = bin edge minimising SSE; 95% CI from
1,000 residual block bootstraps (block = 5 bins). Seed fixed (44).
Outputs: analysis/<RUN>/stage9_supporting/depth_boundary/{depth_bins_200bp.tsv, boundary_estimates.tsv}
"""
import os
import subprocess
import sys
import numpy as np
sys.path.insert(0, os.path.dirname(__file__))
import p44lib as L

RUN = os.environ.get("P44_RUN_DATE", "20261024")
OUT = f"{L.PROJ}/analysis/{RUN}_read_rescue/stage9_supporting/depth_boundary"
os.makedirs(OUT, exist_ok=True)
SAMTOOLS = "/Work_bio/gao/configs/.conda/envs/regular_bioinfo/bin/samtools"
REG = (114580000, 114650000)
NORM = (114400000, 114580000)
BIN = 200
FIT = (114590000, 114640000)
rows = [r for r in L.samples(reported_only=True)]
rng = np.random.default_rng(44)


def depth(bam, start, end, q):
    out = subprocess.run(f"{SAMTOOLS} depth -a -Q {q} -r chr11:{start}-{end} {bam}", shell=True, check=True,
                         capture_output=True, text=True).stdout
    return np.array([int(l.split("\t")[2]) for l in out.splitlines()], dtype=float)


nb = (REG[1] - REG[0] + 1) // BIN
edges = REG[0] + np.arange(nb) * BIN
res = {}
for q in (20, 0):
    for r in rows:
        s = r["sample_id"]
        d = depth(r["bam_path"], *REG, q)
        n = depth(r["bam_path"], *NORM, q)
        med = np.median(n[n > 0])
        bins = d[:nb * BIN].reshape(nb, BIN).mean(axis=1) / med
        res[(q, s)] = (bins, med)
        print(q, s, f"norm_median={med:.1f}", flush=True)

with open(f"{OUT}/depth_bins_200bp.tsv", "w") as fh:
    fh.write("mapq_min\tbin_start\tbin_end\t" + "\t".join(r["sample_id"] for r in rows) + "\n")
    for q in (20, 0):
        for i, e in enumerate(edges):
            fh.write(f"{q}\t{e}\t{e + BIN - 1}\t" + "\t".join(f"{res[(q, r['sample_id'])][0][i]:.3f}" for r in rows) + "\n")


def step_fit(y):
    best = (np.inf, None)
    for k in range(3, len(y) - 3):
        sse = ((y[:k] - y[:k].mean()) ** 2).sum() + ((y[k:] - y[k:].mean()) ** 2).sum()
        if sse < best[0]:
            best = (sse, k)
    k = best[1]
    return k, y[:k].mean(), y[k:].mean()


sel = (edges >= FIT[0]) & (edges + BIN - 1 <= FIT[1])
E = edges[sel]
out = []
groups = {"DLD1": (["S01_DLD1_pop", "S02_DLD1_c1_ctrl", "S03_DLD1_c2_ctrl"], ["S04_DLD1_c3_dox", "S05_DLD1_c4_dox",
                                                                           "S01_DLD1_pop", "S02_DLD1_c1_ctrl", "S03_DLD1_c2_ctrl"]),
          "GM12878": (["S06_GM12878_pop"], ["S07_GM12878_c5", "S08_GM12878_c6", "S09_GM12878_c7", "S10_GM12878_c8"])}
for q in (20, 0):
    for line, (base, tests) in groups.items():
        b = np.mean([res[(q, s)][0] for s in base], axis=0)[sel]
        for s in tests:
            t = res[(q, s)][0][sel]
            ok = (b > 0.05) & (t >= 0)
            y = np.log2((t[ok] + 0.01) / (b[ok] + 0.01))
            x = E[ok]
            k, l1, l2 = step_fit(y)
            fit = np.where(np.arange(len(y)) < k, l1, l2)
            resid = y - fit
            boots = []
            for _ in range(1000):
                starts = rng.integers(0, len(y) - 4, len(y) // 5 + 2)        # full 5-bin blocks only
                idx = np.concatenate([np.arange(i, i + 5) for i in starts])[:len(y)]
                yb = fit + resid[idx]
                boots.append(step_fit(yb)[0])
            lo, hi = np.percentile(boots, [2.5, 97.5]).astype(int)
            out.append((q, line, s, ",".join(base), int(x[k]), int(x[lo]), int(x[min(hi, len(x) - 1)]),
                        round(2 ** l1, 3), round(2 ** l2, 3), round(l2 - l1, 3)))
with open(f"{OUT}/boundary_estimates.tsv", "w") as fh:
    fh.write("mapq_min\tcell_line\tsample\tbaseline\tboundary_bin_start\tci95_lo\tci95_hi\tratio_left\tratio_right\tlog2_step\n")
    for r in out:
        fh.write("\t".join(map(str, r)) + "\n")
for r in out:
    print(r)
print("stage9c OK")
