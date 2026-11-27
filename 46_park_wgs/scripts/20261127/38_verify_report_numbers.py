"""P46 1127 报告数字核验：每个数从源表现读现算，与报告正文逐条比对。
同一事实出现在报告与表两处 ⇒ 必须有一条检查去比对（2026-08-04 proj26 教训）。"""
import csv, re, sys
P = "/home/gao/projects_2026Q4/46_park_wgs"
R = open(f"{P}/custom_research_report_20261127/Dual_Cut_Project_WGS_1127.md").read()
D = f"{P}/custom_research_report_20261127/vector"
def tsv(f): return list(csv.DictReader(open(f"{D}/{f}"), delimiter="\t"))
def col(rows, k, f=float): return {r["sample"]: f(r[k]) for r in rows}

cop = tsv("backbone_copies_single_reference.tsv")
pms = tsv("pMSCV_only_copies.tsv")
jx  = tsv("backbone_junction_crosscheck.tsv")
ctx = tsv("backbone_free_dna_context.tsv")
pick = lambda rows, e: {r["sample"]: float(r["copies_per_diploid"]) for r in rows if r["element"] == e}
ampr, mch, cas9 = pick(cop, "AmpR"), pick(cop, "mCherry"), pick(pms, "Cas9")
nb = {r["sample"]: int(r["junctions_in_backbone"]) for r in jx}
nc = {r["sample"]: int(r["junctions_in_cassette"]) for r in jx}
dup = col(ctx, "dup_pct")

fails = []
def need(text, why):
    if text not in R: fails.append(f"报告里找不到「{text}」 —— {why}")

o = lambda d: sorted(v for k, v in d.items() if k != "Sample_6")
need(f"{o(ampr)[0]:.2f}–{o(ampr)[-1]:.2f} copies (Sample_6: {ampr['Sample_6']:.2f})", "AmpR 范围")
need(f"{min(cas9.values()):.3f}–{max(cas9.values()):.3f} copies", "Cas9 范围")
need(f"{o(mch)[0]:.2f}–{o(mch)[-1]:.2f} copies (Sample_6: {mch['Sample_6']:.2f})", "mCherry 范围")
need(f"{min(nc.values())}–{max(nc.values())} per clone", "盒式接头范围")
need(f"{sum(nc.values())} vector–genome", "接头总数")
lo = min(dup, key=lambda s: dup[s]); hi = max(dup, key=lambda s: dup[s])
need(f"{dup['Sample_5']:.1f}% duplicates with {ampr['Sample_5']:.2f} copies in Sample_5", "低重复率样本")
need(f"{dup['Sample_6']:.1f}% with {ampr['Sample_6']:.2f}", "高重复率样本")
if sum(nb.values()) != 0: fails.append(f"骨架接头应为 0，实为 {sum(nb.values())}")
if "0 in all 10 clones" not in R: fails.append("§7.4 缺「骨架接头 0」断言")
# 分母非空：没有盒式接头就谈不上"骨架 0 条"
if min(nc.values()) <= 0: fails.append("有克隆盒式接头为 0 ⇒ 接头法无检出力，结论不成立")
print("源表实测：")
print(f"  AmpR  其余9样本 {o(ampr)[0]:.3f}–{o(ampr)[-1]:.3f}  Sample_6 {ampr['Sample_6']:.3f}")
print(f"  Cas9  全部10样本 {min(cas9.values()):.3f}–{max(cas9.values()):.3f}")
print(f"  接头  骨架内合计 {sum(nb.values())}  盒式内 {min(nc.values())}–{max(nc.values())} 合计 {sum(nc.values())}")
print(f"  重复率 最低 {lo} {dup[lo]:.1f}%  最高 {hi} {dup[hi]:.1f}%")
if fails:
    print("\nFAIL:"); [print("  -", f) for f in fails]; sys.exit(1)
print("\nALL_NUMBERS_VERIFIED")
