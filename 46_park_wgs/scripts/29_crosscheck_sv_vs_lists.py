"""P46 交付前复核 5 —— DELLY 私有 SV 是否落在脱靶候选位点或载体整合接头附近。

计划 §8.4 要求：落在脱靶候选上的私有事件才有机制解释；否则只能算克隆各自的本底。
私有 SV 定义：仅 1 个样本为非参考基因型（GT 含 1）、其余样本为 0/0。
邻近判据：SV 任一断点与候选位点距离 ≤ 1 kb（脱靶切点）或 ≤ 5 kb（整合接头，接头坐标本身有不确定性）。
输出：private_sv_per_sample.tsv（逐样本计数 + 组间精确置换 p）、private_sv_hits_lists.tsv（命中明细）。
"""
import subprocess, csv, json, sys, os, collections
from itertools import combinations

PROJ = "/home/gao/projects_2026Q4/46_park_wgs"
BCF, OUTD = sys.argv[1], sys.argv[2]
os.makedirs(OUTD, exist_ok=True)
NEAR_OT, NEAR_INT = 1000, 5000

st = list(csv.DictReader(open(f"{PROJ}/scripts/0_sample_table.csv")))
role = {r["safe_name"]: ("control" if r["role"].startswith("control") else "treated") for r in st}
names = subprocess.run(["bcftools", "query", "-l", BCF], capture_output=True, text=True, check=True).stdout.split()
assert set(names) == set(role), (names, list(role))

# 参考清单
ot = collections.defaultdict(list)
for r in csv.DictReader(open(f"{PROJ}/analysis/stage2b_offtarget_evidence_20261123_v2/tables/offtarget_evidence_n10.tsv"), delimiter="\t"):
    if r["is_ontarget"] != "1":
        ot[r["chrom"]].append((int(r["cut_left"]), r["guide"], r["mismatches"]))
integ = collections.defaultdict(list)
for r in csv.DictReader(open(f"{PROJ}/analysis/stage3_vector_v2_20261124/tables/integration_sharing_n10.tsv"), delimiter="\t"):
    integ[r["chrom"]].append((int(r["breakpoint"]), r["n_samples"], r["gene"]))

fmt = "%CHROM\t%POS\t%INFO/SVTYPE\t%INFO/END\t%INFO/CHR2\t%FILTER[\t%GT]\n"
txt = subprocess.run(["bcftools", "query", "-f", fmt, BCF], capture_output=True, text=True, check=True).stdout

def near(chrom, pos, table, win):
    best = None
    for p, *meta in table.get(chrom, []):
        d = abs(p - pos)
        if d <= win and (best is None or d < best[0]): best = (d, p, meta)
    return best

audit = collections.Counter()
priv = collections.Counter(); hits = []
for line in txt.rstrip("\n").split("\n"):
    if not line: continue
    f = line.split("\t")
    chrom, pos, svtype, end, chr2, filt = f[0], int(f[1]), f[2], f[3], f[4], f[5]
    gts = f[6:]
    audit["sv_total"] += 1
    if filt not in ("PASS", "."):
        audit["excl_not_pass"] += 1; continue
    carriers = [s for s, g in zip(names, gts) if "1" in g.replace("|", "/").split("/")]
    homref = sum(1 for g in gts if g.replace("|", "/") == "0/0")
    if len(carriers) != 1:
        audit["excl_not_single_carrier"] += 1; continue
    if homref != len(names) - 1:
        audit["excl_other_samples_not_homref"] += 1; continue
    s = carriers[0]; priv[s] += 1; audit["kept_private_sv"] += 1
    ends = [(chrom, pos)]
    if end not in (".", "") and end.isdigit(): ends.append((chrom, int(end)))
    if chr2 not in (".", ""): ends.append((chr2, int(end) if end.isdigit() else pos))
    for c, p in ends:
        a = near(c, p, ot, NEAR_OT)
        if a: hits.append(dict(sample=s, group=role[s], svtype=svtype, sv_chrom=c, sv_pos=p,
                               list="offtarget_candidate", dist=a[0], site=a[1], detail=",".join(map(str, a[2]))))
        b = near(c, p, integ, NEAR_INT)
        if b: hits.append(dict(sample=s, group=role[s], svtype=svtype, sv_chrom=c, sv_pos=p,
                               list="integration_junction", dist=b[0], site=b[1], detail=",".join(map(str, b[2]))))

with open(f"{OUTD}/private_sv_hits_lists.tsv", "w", newline="") as fh:
    cols = ["sample", "group", "svtype", "sv_chrom", "sv_pos", "list", "dist", "site", "detail"]
    w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t"); w.writeheader(); w.writerows(hits)

def exact_p(a, b):
    allv = a + b; n = len(allv); k = len(b)
    obs = abs(sum(a) / len(a) - sum(b) / len(b)); c = t = 0
    for idx in combinations(range(n), k):
        bb = [allv[i] for i in idx]; aa = [allv[i] for i in range(n) if i not in idx]; t += 1
        if abs(sum(aa) / len(aa) - sum(bb) / len(bb)) >= obs - 1e-12: c += 1
    return c / t

per = [dict(sample=s, group=role[s], private_sv=priv[s],
            hits_offtarget=sum(1 for h in hits if h["sample"] == s and h["list"] == "offtarget_candidate"),
            hits_integration=sum(1 for h in hits if h["sample"] == s and h["list"] == "integration_junction"))
       for s in names]
with open(f"{OUTD}/private_sv_per_sample.tsv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(per[0]), delimiter="\t"); w.writeheader(); w.writerows(per)

t = [d["private_sv"] for d in per if d["group"] == "treated"]
c = [d["private_sv"] for d in per if d["group"] == "control"]
summary = dict(audit=dict(audit), per_sample=per,
               treated=t, control=c, p_exact_two_sided=round(exact_p(t, c), 4),
               min_attainable_p=round(1 / 120, 4),
               near_window=dict(offtarget_bp=NEAR_OT, integration_bp=NEAR_INT),
               constraint="SV 仅 DELLY 单一方法（Manta 在本数据集不可行）⇒ 断点级结论按『候选』表述；"
                          "且编辑组与对照组核型不同源（chr12p），组间差异只能表述为『与编辑共分离』。")
json.dump(summary, open(f"{OUTD}/private_sv_summary.json", "w"), indent=1, ensure_ascii=False)
print(json.dumps(summary, indent=1, ensure_ascii=False))
