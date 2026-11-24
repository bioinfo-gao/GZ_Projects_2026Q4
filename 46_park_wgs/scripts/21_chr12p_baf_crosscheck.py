"""P46 — chr12p 拷贝数差异的独立互证：杂合位点 B 等位比例（BAF）。

为什么需要：Stage 4c 的臂级拷贝估计（深度法）给出 chr12p 处理组 ~4 拷贝 / 对照组 ~2 拷贝。
深度法与 BAF 法的失效模式完全不同（深度受文库/GC/重复区影响，BAF 只看同一位点内两个等位的读段比），
⇒ 两者一致才足以写进报告（feedback_self_audit_needs_independent_cross_check）。

预期：2 拷贝杂合 ⇒ BAF ≈ 0.50；4 拷贝（3:1）⇒ 出现 BAF ≈ 0.25 的一群；12q 两组都应 ≈ 0.50（内部对照）。
输入复用身份矩阵那一轮的 10 样本联合 pileup（gnomAD AF 0.3–0.7 常见位点，已带 AD）。
"""
import subprocess, statistics, collections, csv, json, os, sys

PROJ = "/home/gao/projects_2026Q4/46_park_wgs"
V = f"{PROJ}/analysis/stage0_identity_qc_20261123_v2/pileup/all.vcf.gz"
OUT = f"{PROJ}/analysis/stage4c_cnv_20261123/baf_check"
CEN = 35_500_000          # chr12 着丝粒近似位置（GRCh38）
MIN_DP = 20
os.makedirs(OUT, exist_ok=True)

names = subprocess.run(["bcftools", "query", "-l", V], capture_output=True, text=True, check=True).stdout.split()
assert len(names) == 10, names
st = {r["safe_name"]: ("control" if r["role"].startswith("control") else "treated")
      for r in csv.DictReader(open(f"{PROJ}/scripts/0_sample_table.csv"))}

def load(region):
    txt = subprocess.run(["bcftools", "query", "-f", "%CHROM\t%POS[\t%GT:%AD]\n", "-r", region, V],
                         capture_output=True, text=True, check=True).stdout
    return [ln.split("\t") for ln in txt.rstrip("\n").split("\n") if ln]

groups = {"chr12p": [], "chr12q": [], "chr7(对照臂)": []}
for f in load("chr12"):
    groups["chr12p" if int(f[1]) < CEN else "chr12q"].append(f[2:])
for f in load("chr7"):
    groups["chr7(对照臂)"].append(f[2:])
print("位点数:", {k: len(v) for k, v in groups.items()})

rows = []
for i, s in enumerate(names):
    for arm, cells in groups.items():
        baf = []
        for c in cells:
            gt, _, ad = c[i].partition(":")
            if gt.replace("|", "/") != "0/1" or ad in ("", "."):
                continue
            v = [int(x) for x in ad.split(",") if x != "."]
            if len(v) < 2 or sum(v) < MIN_DP:
                continue
            baf.append(min(v) / sum(v))
        if not baf:
            continue
        low = [b for b in baf if b < 0.42]
        rows.append(dict(sample=s, group=st[s], arm=arm, n_het=len(baf),
                         median_baf=round(statistics.median(baf), 3),
                         frac_below_0_40=round(sum(1 for b in baf if b < 0.40) / len(baf), 3),
                         median_of_low_mode=round(statistics.median(low), 3) if low else None))
with open(f"{OUT}/chr12p_baf_crosscheck.tsv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t"); w.writeheader(); w.writerows(rows)

print(f"{'sample':<10}{'group':<9}{'arm':<14}{'n_het':>7}{'medBAF':>8}{'frac<0.40':>10}{'lowMode':>9}")
for r in rows:
    print(f"{r['sample']:<10}{r['group']:<9}{r['arm']:<14}{r['n_het']:>7}{r['median_baf']:>8}"
          f"{r['frac_below_0_40']:>10}{str(r['median_of_low_mode']):>9}")

# 判据：12p 上处理组 frac<0.40 应显著高于对照组；12q 与 chr7 两组应无差别
def fr(arm, grp): return [r["frac_below_0_40"] for r in rows if r["arm"] == arm and r["group"] == grp]
verdict = {}
for arm in groups:
    t, c = fr(arm, "treated"), fr(arm, "control")
    verdict[arm] = dict(treated_min=min(t), treated_max=max(t), control_min=min(c), control_max=max(c),
                        separated=bool(min(t) > max(c)))
print(json.dumps(verdict, indent=1, ensure_ascii=False))
json.dump(verdict, open(f"{OUT}/verdict.json", "w"), indent=1, ensure_ascii=False)
print("[concordance] 深度法(4c)说 12p 处理组 ~4 拷贝；BAF 法若 separated=True 则两法一致")
