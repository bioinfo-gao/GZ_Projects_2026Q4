"""P46 v3b 汇总（2026-11-27）：pMSCV 单载体索引的元件拷贝数。

目的：让 Fig 6 两个 panel 的方法对称（都用单载体索引），并给出 v2(双载体) 与 v3b(单载体) 的
一致性检查 —— Cas9/EGFP 在 pLPC 里无同源序列（blastn 验证），两种口径应当给出一致结果；
若不一致，说明某一侧另有同源竞争，必须查清后才能用。
"""
import csv, os, sys

P = "/home/gao/projects_2026Q4/46_park_wgs"
V2 = f"{P}/analysis/stage3_vector_v2_20261124"
OUT = f"{P}/analysis/stage3_pMSCVonly_v3_20261127"
QC = "/srv/external_data/gao_projects/46_park_wgs_bam/stage0_20261121/qc"
G = 3215721559.0

NAMES = [r["safe_name"] for r in csv.DictReader(open(f"{P}/scripts/0_sample_table.csv"))]
ROLE = {r["safe_name"]: ("control" if r["role"].startswith("control") else "treated")
        for r in csv.DictReader(open(f"{P}/scripts/0_sample_table.csv"))}
ELEM = [("WHOLE_VECTOR", 1, 11191), ("Cas9", 1476, 5576), ("EGFP", 5640, 6356),
        ("PGK_promoter", 6367, 6866), ("shared_backbone", 8407, 10390)]

def genome_depth(S):
    v = {}
    for ln in open(f"{QC}/{S}.samtools_stats.txt"):
        p = ln.rstrip("\n").split("\t")
        if p[0] == "SN":
            v[p[1]] = float(p[2])
    return (v["bases mapped (cigar):"] - v["bases duplicated:"]) / G

def depth_map(path, want=None):
    d = {}
    for ln in open(path):
        f = ln.split()
        if want is None or f[0] == want:
            d[int(f[1])] = int(f[2])
    return d

rows, check = [], []
for S in NAMES:
    f = f"{OUT}/tables/{S}.pMSCVonly_depth.txt"
    if not os.path.exists(f):
        sys.exit(f"FATAL 缺少 {f}")
    d = depth_map(f)
    gd = genome_depth(S)
    d2 = depth_map(f"{V2}/tables/{S}.vec_depth.txt", want="pMSCV_Cas9_EGFP")
    for nm, a, b in ELEM:
        n = b - a + 1
        mean = sum(d.get(p, 0) for p in range(a, b + 1)) / float(n)
        cov = 100.0 * sum(1 for p in range(a, b + 1) if d.get(p, 0) > 0) / n
        rows.append(dict(sample=S, group=ROLE[S], element=nm, start=a, end=b, length_bp=n,
                         mean_depth=round(mean, 2), covered_pct=round(cov, 1),
                         genome_depth=round(gd, 2), copies_per_diploid=round(2 * mean / gd, 3)))
        if nm in ("Cas9", "EGFP"):
            m2 = sum(d2.get(p, 0) for p in range(a, b + 1)) / float(n)
            check.append(dict(sample=S, element=nm, v2_dual_copies=round(2 * m2 / gd, 3),
                              v3b_single_copies=round(2 * mean / gd, 3)))

with open(f"{OUT}/tables/pMSCVonly_copies.tsv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
    w.writeheader(); w.writerows(rows)
with open(f"{OUT}/tables/v2_vs_v3b_consistency.tsv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(check[0]), delimiter="\t")
    w.writeheader(); w.writerows(check)

# 一致性断言：Cas9/EGFP 无同源竞争 ⇒ 两种口径差异必须很小
bad = [c for c in check
       if abs(c["v2_dual_copies"] - c["v3b_single_copies"]) > 0.05 + 0.10 * max(c["v2_dual_copies"], c["v3b_single_copies"])]
print("V3B_SUMMARY_DONE")
for c in check:
    print(f"  {c['sample']:9s} {c['element']:5s} v2双载体={c['v2_dual_copies']:>8.3f}  v3b单载体={c['v3b_single_copies']:>8.3f}")
if bad:
    sys.exit("FATAL 两种口径不一致，说明另有同源竞争，必须查清：%s" % bad)
print("CONSISTENCY_OK  Cas9/EGFP 在两种口径下一致（10% + 0.05 拷贝容差）")
