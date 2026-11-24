"""10×10 胚系一致率 + hom-alt 位点 REF 读段占比（外源污染指标）。
输入：bcftools query -H '%CHROM\t%POS[\t%GT:%AD]'；位点仅计两样本 DP≥10 都有 GT 的。
输出：identity_matrix.tsv（一致率）、identity_n_sites.tsv（分母）、contamination.tsv、filter_impact.tsv。
"""
import sys, itertools, csv
inp, outd = sys.argv[1], sys.argv[2]
MIN_DP = 10
rows = []
with open(inp) as f:
    hdr = f.readline().rstrip("\n").split("\t")
    names = [h.split("]")[1].split(":")[0].replace(".markdup.bam", "").split("/")[-1] for h in hdr[2:]]
    for line in f:
        rows.append(line.rstrip("\n").split("\t")[2:])
n = len(names)
assert n == 10, names

def parse(x):
    gt, _, ad = x.partition(":")
    if "." in gt: return None, 0, None
    a = [int(v) for v in ad.split(",") if v != "."] if ad and ad != "." else []
    dp = sum(a)
    g = gt.replace("|", "/")
    return (g if dp >= MIN_DP else None), dp, a

audit = {s: {"sites_total": 0, "no_call": 0, "dp_lt_10": 0, "used": 0} for s in names}
P = []
for r in rows:
    pr = []
    for s, x in zip(names, r):
        audit[s]["sites_total"] += 1
        g, dp, a = parse(x)
        if "." in x.split(":")[0]: audit[s]["no_call"] += 1
        elif dp < MIN_DP: audit[s]["dp_lt_10"] += 1
        else: audit[s]["used"] += 1
        pr.append((g, a))
    P.append(pr)

conc = [[0.0]*n for _ in range(n)]; nsite = [[0]*n for _ in range(n)]
for i, j in itertools.product(range(n), range(n)):
    m = t = 0
    for pr in P:
        gi, gj = pr[i][0], pr[j][0]
        if gi is None or gj is None: continue
        t += 1; m += (gi == gj)
    conc[i][j] = m / t if t else float("nan"); nsite[i][j] = t

def wmat(fn, M, fmt):
    with open(f"{outd}/{fn}", "w") as o:
        o.write("sample\t" + "\t".join(names) + "\n")
        for i, s in enumerate(names): o.write(s + "\t" + "\t".join(fmt(v) for v in M[i]) + "\n")
wmat("identity_matrix.tsv", conc, lambda v: f"{v:.4f}")
wmat("identity_n_sites.tsv", nsite, str)

with open(f"{outd}/contamination.tsv", "w") as o:
    o.write("sample\thom_alt_sites\tref_reads\ttotal_reads\tref_frac_at_homalt\thet_sites\thet_alt_frac_median\n")
    for i, s in enumerate(names):
        hs = rr = tt = 0; het = []
        for pr in P:
            g, a = pr[i]
            if g is None or len(a) < 2: continue
            if g == "1/1": hs += 1; rr += a[0]; tt += sum(a)
            elif g == "0/1": het.append(a[1] / sum(a))
        het.sort()
        med = het[len(het)//2] if het else float("nan")
        o.write(f"{s}\t{hs}\t{rr}\t{tt}\t{rr/tt if tt else float('nan'):.5f}\t{len(het)}\t{med:.3f}\n")

with open(f"{outd}/filter_impact.tsv", "w") as o:
    o.write("sample\tsites_total\tno_call\tdp_lt_10\tused\n")
    for s in names: o.write(s + "\t" + "\t".join(str(audit[s][k]) for k in ("sites_total","no_call","dp_lt_10","used")) + "\n")

off = [conc[i][j] for i in range(n) for j in range(n) if i != j]
print(f"[identity] off-diagonal concordance min={min(off):.4f} max={max(off):.4f}; min shared sites={min(nsite[i][j] for i in range(n) for j in range(n) if i!=j)}")
