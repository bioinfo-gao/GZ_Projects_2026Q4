"""P46 Stage 3 v3 汇总（2026-11-27）。

两项彼此独立、失效模式不同的测量：
  主测 —— 候选读段回贴【单载体】pLPC 索引（无同源竞争者），`-Q 20` 不再把骨架读段压成 MAPQ 0。
  互证 —— 在 pLPC/pMSCV 的**分歧位点**上数等位：骨架 DNA 到底来自哪个质粒。
          用【双载体】BAM 且**不设 MAPQ 下限**（正是 v2 丢掉的那批读段）。

两者只要有一项说"骨架不存在"，就不能写成"存在"；两者都说存在才下结论。
拷贝数只能回答"在不在、多少"，**整合 vs 游离污染另需接头证据**，本脚本不下该结论。
"""
import csv, os, sys, subprocess, collections, json

P = "/home/gao/projects_2026Q4/46_park_wgs"
V2 = f"{P}/analysis/stage3_vector_v2_20261124"
OUT = f"{P}/analysis/stage3_backbone_v3_20261127"
REFS = f"{P}/analysis/stage3_refs_20261121"
QC = "/srv/external_data/gao_projects/46_park_wgs_bam/stage0_20261121/qc"
G = 3215721559.0
os.makedirs(f"{OUT}/tables", exist_ok=True)

NAMES = [r["safe_name"] for r in csv.DictReader(open(f"{P}/scripts/0_sample_table.csv"))]
ROLE = {r["safe_name"]: ("control" if r["role"].startswith("control") else "treated")
        for r in csv.DictReader(open(f"{P}/scripts/0_sample_table.csv"))}

def genome_depth(S):
    v = {}
    for ln in open(f"{QC}/{S}.samtools_stats.txt"):
        p = ln.rstrip("\n").split("\t")
        if p[0] == "SN":
            v[p[1]] = float(p[2])
    return (v["bases mapped (cigar):"] - v["bases duplicated:"]) / G

# pLPC 元件（与 v2 同一来源，只取需要解释的几个）
ELEM = [("WHOLE_VECTOR", 1, 7985), ("ori", 234, 822), ("AmpR", 982, 1842),
        ("PuroR", 3346, 3945), ("mCherry", 4956, 5663), ("LTR_5p", 1978, 2447),
        ("LTR_3p", 7449, 7918)]
# 53BP1 区段与人源同源，v2 已 mask；此处同样排除，口径一致
MASK = set(range(5550, 7261))

# ── 主测：单载体索引深度 ────────────────────────────────────────────
rows = []
for S in NAMES:
    f = f"{OUT}/tables/{S}.pLPConly_depth.txt"
    if not os.path.exists(f):
        sys.exit(f"FATAL 缺少 {f}（映射未完成？）")
    d = {}
    for ln in open(f):
        _, p, x = ln.split()
        d[int(p)] = int(x)
    gd = genome_depth(S)
    for nm, a, b in ELEM:
        pos = [p for p in range(a, b + 1) if p not in MASK]
        mean = sum(d.get(p, 0) for p in pos) / len(pos)
        cov = 100.0 * sum(1 for p in pos if d.get(p, 0) > 0) / len(pos)
        rows.append(dict(sample=S, group=ROLE[S], element=nm, start=a, end=b,
                         unmasked_bp=len(pos), mean_depth=round(mean, 2),
                         covered_pct=round(cov, 1), genome_depth=round(gd, 2),
                         copies_per_diploid=round(2 * mean / gd, 3)))
with open(f"{OUT}/tables/backbone_v3_copies.tsv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t")
    w.writeheader(); w.writerows(rows)

# ── 互证：整合接头的载体端落位（与深度法失效模式完全不同）──────────────
# 分歧位点法已废弃：pLPC 1-1973 与 pMSCV 8407-10390 的 99.1% 差异里只有 7 个替换、
# 落在 ori/AmpR 内的仅 1 个 ⇒ 检出力不足，不能用"0 条"下任何结论（见 research record）。
# 改用接头证据：逆转录病毒载体只有 LTR 之间的盒式序列会整合，质粒骨架不会。
# 若骨架确实整合，必然出现骨架坐标上的 vector–genome 接头。
J = list(csv.DictReader(open(f"{V2}/tables/integration_junctions_n10.tsv"), delimiter="\t"))
BB = (1, 1973)                       # 与 pMSCV 共享的 ori+AmpR 骨架
CASSETTE = (1978, 7918)              # 5'LTR 到 3'LTR，即会整合的那一段
jx = []
for S in NAMES:
    js = [r for r in J if r["sample"] == S and "pLPC" in r["vector"]]
    nb = sum(1 for r in js if BB[0] <= int(r["vec_pos_min"]) <= BB[1])
    nc = sum(1 for r in js if CASSETTE[0] <= int(r["vec_pos_min"]) <= CASSETTE[1])
    jx.append(dict(sample=S, group=ROLE[S], junctions_in_backbone=nb,
                   junctions_in_cassette=nc, total_pLPC_junctions=len(js)))
with open(f"{OUT}/tables/backbone_v3_junction_crosscheck.tsv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(jx[0]), delimiter="\t")
    w.writeheader(); w.writerows(jx)
# 分母非空断言：盒式序列必须有接头，否则"骨架 0 条"同样可能是没有检出力
assert sum(x["junctions_in_cassette"] for x in jx) > 0, "盒式序列接头为 0 ⇒ 接头法本身无检出力，不得据此下结论"

# 游离 DNA vs 整合：骨架拷贝数与文库重复率的关系
def dup_rate(S):
    fl = open(f"{QC}/{S}.flagstat.txt").read().split("\n")
    pr = int([l for l in fl if l.endswith(" primary")][0].split()[0])
    du = int([l for l in fl if "primary duplicates" in l][0].split()[0])
    return 100.0 * du / pr
ampr = {r["sample"]: r["copies_per_diploid"] for r in rows if r["element"] == "AmpR"}
mch = {r["sample"]: r["copies_per_diploid"] for r in rows if r["element"] == "mCherry"}
ctx = [dict(sample=S, group=ROLE[S], dup_pct=round(dup_rate(S), 1),
            AmpR_copies=ampr[S], mCherry_copies=mch[S]) for S in NAMES]
with open(f"{OUT}/tables/backbone_v3_context.tsv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(ctx[0]), delimiter="\t")
    w.writeheader(); w.writerows(ctx)

# ── 审计：v2 口径 vs v3 口径，同一区段 ────────────────────────────────
aud = []
for S in NAMES:
    d2 = {}
    for ln in open(f"{V2}/tables/{S}.vec_depth.txt"):
        v, p, x = ln.split()
        if v == "pLPC_mCherry_BP1_2":
            d2[int(p)] = int(x)
    v2_ampr = sum(d2.get(p, 0) for p in range(982, 1843)) / 861.0
    r = [x for x in rows if x["sample"] == S and x["element"] == "AmpR"][0]
    n_reads = subprocess.run(["samtools", "view", "-c", f"{V2}/bam/{S}.vec.bam",
                              "pLPC_mCherry_BP1_2:982-1842"], capture_output=True, text=True).stdout.strip()
    n_q20 = subprocess.run(["samtools", "view", "-c", "-q", "20", f"{V2}/bam/{S}.vec.bam",
                            "pLPC_mCherry_BP1_2:982-1842"], capture_output=True, text=True).stdout.strip()
    aud.append(dict(sample=S, v2_AmpR_mean_depth=round(v2_ampr, 2),
                    v2_reads_in_AmpR=n_reads, v2_reads_mapq_ge20=n_q20,
                    v3_AmpR_mean_depth=r["mean_depth"], v3_AmpR_covered_pct=r["covered_pct"],
                    v3_AmpR_copies=r["copies_per_diploid"]))
with open(f"{OUT}/tables/backbone_v2_vs_v3_audit.tsv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(aud[0]), delimiter="\t")
    w.writeheader(); w.writerows(aud)

json.dump(dict(backbone_region_pLPC=list(BB), cassette_region_pLPC=list(CASSETTE),
               junctions_in_backbone_total=sum(x["junctions_in_backbone"] for x in jx),
               junctions_in_cassette_total=sum(x["junctions_in_cassette"] for x in jx)),
          open(f"{OUT}/tables/backbone_v3_meta.json", "w"), indent=1)
print("V3_SUMMARY_DONE  骨架接头合计 =", sum(x["junctions_in_backbone"] for x in jx),
      " 盒式接头合计 =", sum(x["junctions_in_cassette"] for x in jx))
for a in aud:
    print(f"  {a['sample']:9s} v2深度={a['v2_AmpR_mean_depth']:6.2f}(读段{a['v2_reads_in_AmpR']:>5s},MAPQ≥20 {a['v2_reads_mapq_ge20']:>4s})"
          f"   v3深度={a['v3_AmpR_mean_depth']:7.2f} 覆盖{a['v3_AmpR_covered_pct']:5.1f}% 拷贝={a['v3_AmpR_copies']}")
for x in jx:
    print(f"  {x['sample']:9s} 接头: 骨架内={x['junctions_in_backbone']:>3d}  盒式内={x['junctions_in_cassette']:>3d}  pLPC合计={x['total_pLPC_junctions']:>3d}")
for c in ctx:
    print(f"  {c['sample']:9s} dup={c['dup_pct']:5.1f}%  AmpR={c['AmpR_copies']:>7.3f}  mCherry={c['mCherry_copies']:>7.3f}")
