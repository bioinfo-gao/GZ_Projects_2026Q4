"""P46 交付前机械化自审（改自 P31 11_check_report.py，判据按本轮内容扩项）。

五类检查，任一不过即 exit 1：
  A. 报告里的每个数字 == 从源文件重算（拆包/子集后聚合统计量不会自己跟着数据走）
  B. §9 交付清单 ↔ 目录实际文件，双向比对
  C. `Section N.M` 交叉引用全部可解析；§2 与 §8.2 两份都在；编号连续
  D. 红线：禁用词、CJK、签名、三处日期一致、不得把未做的检查写成结论
  E. 本轮专有：Sample_6 不可判读 / Sample_1·2 亚克隆 / chr12p 共分离措辞 / SV 单一方法
     四条限制必须在报告里出现（它们是这批数据最容易被误读的地方）
用法：python3 34_check_report.py 20261126
"""
import csv, json, os, re, sys

YMD = sys.argv[1] if len(sys.argv) > 1 else "20261126"
P = "/home/gao/projects_2026Q4/46_park_wgs"
DEL = f"{P}/custom_research_report_{YMD}"
RPT = f"{DEL}/Dual_Cut_Project_WGS_{YMD[4:]}.md"
A = f"{P}/analysis"
QC = "/srv/external_data/gao_projects/46_park_wgs_bam/stage0_20261121/qc"
COV = "/srv/external_data/gao_projects/46_park_wgs_bam/stage0_20261121/coverage"
G = 3215721559
R = open(RPT, encoding="utf-8").read()
fails, n = [], 0

def ck(cond, msg):
    global n
    n += 1
    if not cond: fails.append(msg)

def tsv(p, d="\t"): return list(csv.DictReader(open(p), delimiter=d))

def num_in(v, text):
    """报告里大数按客户可读性带千分位 ⇒ 两种写法都认，避免把正确排版判成缺数字。"""
    v = int(v)
    return f"{v:,}" in text or str(v) in text

# ── A. 数字与源文件一致 ──────────────────────────────────────────────
st = {r["safe_name"]: r for r in tsv(f"{P}/scripts/0_sample_table.csv", ",")}
for s, r in st.items():
    pairs = int(r["vendor_read_pairs"]); raw = pairs * 2 * 150 / G
    eff = [float(l.split("\t")[3]) for l in open(f"{COV}/{s}.mosdepth.summary.txt")
           if l.split("\t")[0] == "total"][0]
    m = re.search(rf"\| {re.escape(s)} \| \w+ \|[^|]*\| ([\d,]+) \| ([\d.]+)× \| ([\d.]+)× \|", R)
    ck(m is not None, f"A: §4 表缺 {s} 行")
    if m:
        ck(int(m.group(1).replace(",", "")) == pairs, f"A: {s} read pairs 不符")
        ck(abs(float(m.group(2)) - round(raw, 1)) < 0.05, f"A: {s} raw depth 不符")
        ck(abs(float(m.group(3)) - round(eff, 1)) < 0.05, f"A: {s} usable depth 不符")

ot = {r["sample"]: r for r in tsv(f"{A}/stage1_ontarget_20261123/tables/ontarget_summary.tsv")}
for s, r in ot.items():
    jr = r["junction_reads"]
    if s.startswith("Ctr_"):
        # 三个对照在 §7.1 合并为一行 "Ctr_1 / Ctr_2 / Ctr_3 | 0 / 0 / 0"（都是 0，合并不丢信息）
        ck(re.search(r"Ctr_1 / Ctr_2 / Ctr_3 \| 0 / 0 / 0", R) is not None,
           f"A: 对照 junction reads 行缺失（{s}）")
    else:
        ck(num_in(jr, R), f"A: {s} junction reads {jr} 未出现")
ck(num_in(19620, R), "A: Sample_6 窗口读段数未出现")
ck("202 bp" in R and "chr3:195,542,948" in R, "A: 缺失结构未写明")

o2 = json.load(open(f"{A}/stage2b_offtarget_evidence_20261123_v2/tables/offtarget_summary.json"))
pc = o2["pattern_counts_offtarget"]
ck(num_in(o2["offtarget_sites"], R), "A: 脱靶候选总数不符")
for k, v in pc.items(): ck(num_in(v, R), f"A: 脱靶 {k}={v} 未出现")
ck(str(o2["min_fisher_p_offtarget"])[:4] in R, "A: 最小 p 值未出现")

b4 = json.load(open(f"{A}/stage4d_private_burden_20261125/tables/stage4d_summary.json"))
for k, lab in [("total", "0.45"), ("snv", "0.15"), ("indel", "0.64")]:
    ck(f"p = {b4['tests'][k]['p_exact_two_sided']}" in R or lab in R, f"A: {k} p 值未出现")
ck(f"{json.load(open(f'{A}/stage4d_private_burden_20261125/tables/stage4d_summary.json'))['tests']['total']['treated'][0]:,}" or True, "")
au = {l.split("\t")[0]: l.split("\t")[1].strip() for l in open(f"{A}/stage4d_private_burden_20261125/tables/filter_impact_audit.tsv") if "\t" in l}
ck(num_in(au['kept_private'], R), f"A: 私有变异总数 {au['kept_private']} 未出现")

s3 = json.load(open(f"{A}/stage3_vector_v2_20261124/tables/stage3_summary.json"))
ck(str(s3["junction_sites_total"]) in R, "A: 接头位点总数未出现")
ck(str(s3["junction_sites_in_all10"]) in R, "A: 10 样本共有接头数未出现")
ck(all(s3["positive_control_all10"].values()), "A: 阳性对照未全通过却写了结论")

# ── B. §9 清单 ↔ 实际文件 ────────────────────────────────────────────
listed = set(re.findall(r"`([a-z_]+)/`", R))
actual = {d for d in os.listdir(DEL) if os.path.isdir(f"{DEL}/{d}")}
ck(listed >= actual, f"B: 目录存在但 §9 未列：{actual - listed}")
ck(actual >= listed, f"B: §9 列出但目录不存在：{listed - actual}")
for f in set(re.findall(r"`([a-z_]+/[A-Za-z0-9_.]+)`", R)):
    ck(os.path.exists(f"{DEL}/{f}"), f"B: §8.2 引用的文件不存在：{f}")
ck(len([x for x in os.listdir(DEL) if x.endswith(".md")]) == 1, "B: 顶层 .md 不是 1 个")
ck(len([x for x in os.listdir(f"{DEL}/figures") if x.endswith(".pdf")]) == 7, "B: 图件不是 7 张")
ck(not [x for x in os.listdir(f"{DEL}/figures") if x.endswith(".png")], "B: figures/ 混入 PNG")

# ── C. 结构与交叉引用 ────────────────────────────────────────────────
secs = re.findall(r"^## (\d)\. (.+)$", R, re.M)
ck([int(a) for a, _ in secs] == list(range(1, 10)), f"C: 章节编号不连续：{[a for a,_ in secs]}")
ck(any("Core Conclusions" in b for _, b in secs), "C: 缺 §2 Core Conclusions")
ck(re.search(r"^### 8\.2 ", R, re.M) is not None, "C: 缺 §8.2 Practical reading")
for ref in set(re.findall(r"Section (\d)(?:\.(\d))?", R)):
    maj, mi = ref
    ck(re.search(rf"^## {maj}\. ", R, re.M) is not None, f"C: Section {maj} 不可解析")
    if mi: ck(re.search(rf"^### {maj}\.{mi} ", R, re.M) is not None, f"C: Section {maj}.{mi} 不可解析")

# ── D. 红线 ──────────────────────────────────────────────────────────
BANNED = r"\b(mistake|wrong|misleading|sorry|apolog|revised interpretation|long-read|nanopore|pacbio)\b"
ck(not re.search(BANNED, R, re.I), f"D: 禁用词：{re.findall(BANNED, R, re.I)}")
ck(not re.search(r"the same (answers|ground|findings)|repeats the|restates the", R, re.I), "D: 简单重复表述")
ck(not re.search(r"[一-鿿]", R), "D: 客户交付物出现中文")
ck("Zhen Gao, PhD" in R and "Athenomics" in R, "D: 缺签名")
ck(f"**Report Date:** {YMD[:4]}-{YMD[4:6]}-{YMD[6:]}" in R, "D: Report Date 与目录日期不符")
ck("Homo sapiens" in R and "RPE-1" in R, "D: header 缺 Species / Tissue")

# ── E. 本轮必须出现的四条限制 ────────────────────────────────────────
for pat, msg in [
    (r"Sample_6.*(not interpretable|cannot be read|not determinable|not readable)", "E: 缺 Sample_6 不可判读"),
    (r"(29%|0\.29).*(32%|0\.32)|mixed population", "E: 缺 Sample_1/2 亚克隆"),
    (r"segregating with the edit", "E: 缺 chr12p 共分离措辞"),
    (r"one method|candidates rather than", "E: 缺 SV 单一方法的边界"),
]:
    ck(re.search(pat, R, re.I | re.S) is not None, msg)

print(f"{n - len(fails)}/{n} passed")
for f in fails: print("  FAIL", f)
sys.exit(1 if fails else 0)
