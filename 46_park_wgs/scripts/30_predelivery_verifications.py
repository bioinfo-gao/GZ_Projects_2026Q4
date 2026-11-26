"""P46 交付前复核 ①③④⑤（2026-11-26）—— 全部写进新目录，不覆盖任何既有产物。

① 7 个 treated_only 脱靶位点：逐样本导出切点 ±10 bp 内的 indel / clip 读段（含 CIGAR、MAPQ、
   read name、是否 supplementary），判断是否"同一个 indel 在多样本重复出现"（真脱靶的特征）
   还是各样本各自一条（本底特征）。
⑤ Sample_4 的 chr19 INV（唯一一条私有 SV 落在脱靶候选附近）：导出两个断点与候选切点附近的读段。
③ 载体拷贝数按候选读段总数归一：现值与文库重复率同序 ⇒ 未归一的组间比较不可用。
④ Stage 1 定标表：原表把 Sample_6 的"不可判读"计为 agrees=yes（显示 10/10），
   改写为 9/10 一致 + 1 未定。**原表不动**，新表另存并在 README 写明关系。
"""
import csv, json, os, re, subprocess, collections

PROJ = "/home/gao/projects_2026Q4/46_park_wgs"
BAM = "/srv/external_data/gao_projects/46_park_wgs_bam/stage0_20261121/bam"
QC = "/srv/external_data/gao_projects/46_park_wgs_bam/stage0_20261121/qc"
S2B = f"{PROJ}/analysis/stage2b_offtarget_evidence_20261123_v2/tables"
S3 = f"{PROJ}/analysis/stage3_vector_v2_20261124/tables"
S1 = f"{PROJ}/analysis/stage1_ontarget_20261123/tables"
OUT = f"{PROJ}/analysis/stage5_verify_20261126"
FLANK, MIN_CLIP = 10, 20
os.makedirs(f"{OUT}/tables", exist_ok=True)

names = [r["safe_name"] for r in csv.DictReader(open(f"{PROJ}/scripts/0_sample_table.csv"))]
role = {r["safe_name"]: ("control" if r["role"].startswith("control") else "treated")
        for r in csv.DictReader(open(f"{PROJ}/scripts/0_sample_table.csv"))}

def cig(c): return [(int(n), op) for n, op in re.findall(r"(\d+)([MIDNSHP=X])", c)]
def reflen(c): return sum(n for n, op in cig(c) if op in "MDN=X")

# ───────── ① 7 个 treated_only 位点的读段级证据 ─────────
sites = [r for r in csv.DictReader(open(f"{S2B}/offtarget_evidence_n10.tsv"), delimiter="\t")
         if r["is_ontarget"] != "1" and r["pattern"] == "treated_only"]
rows = []
for s in sites:
    chrom, cut = s["chrom"], int(s["cut_left"])
    for smp in names:
        out = subprocess.run(["samtools", "view", f"{BAM}/{smp}.markdup.bam",
                              f"{chrom}:{cut-200}-{cut+200}"], capture_output=True, text=True).stdout
        for ln in out.rstrip("\n").split("\n"):
            if not ln: continue
            f = ln.split("\t")
            flag, pos, mapq, cigar, seq = int(f[1]), int(f[3]), int(f[4]), f[5], f[9]
            if flag & 0x400 or cigar == "*": continue
            end = pos + reflen(cigar) - 1
            if end < cut - FLANK or pos > cut + FLANK + 1: continue
            ops = cig(cigar)
            ev, detail = [], []
            q = pos
            for n, op in ops:
                if op in "ID" and abs(q - cut) <= FLANK:
                    ev.append("indel"); detail.append(f"{op}{n}@{q}")
                if op in "MDN=X": q += n
            adapter = ("AGATCGGAAG" in seq) or ("CTTCCGATCT" in seq)
            if ops[0][1] == "S" and ops[0][0] >= MIN_CLIP and abs(pos - 1 - cut) <= FLANK:
                ev.append("clipL"); detail.append(f"S{ops[0][0]}L")
            if ops[-1][1] == "S" and ops[-1][0] >= MIN_CLIP and abs(end - cut) <= FLANK:
                ev.append("clipR"); detail.append(f"S{ops[-1][0]}R")
            if not ev: continue
            rows.append(dict(guide=s["guide"], chrom=chrom, cut_left=cut, sample=smp, group=role[smp],
                             evidence="+".join(sorted(set(ev))), detail=";".join(detail),
                             cigar=cigar, mapq=mapq, pos=pos,
                             flag_secondary=int(bool(flag & 0x100)), flag_supplementary=int(bool(flag & 0x800)),
                             adapter_seq_in_read=int(adapter), read=f[0]))
with open(f"{OUT}/tables/treated_only_sites_reads.tsv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0]) if rows else ["guide"], delimiter="\t")
    w.writeheader(); w.writerows(rows)

# 同一 indel 是否在多样本重复（真脱靶特征）
rep = collections.defaultdict(set)
for r in rows:
    for d in r["detail"].split(";"):
        if d[0] in "ID": rep[(r["chrom"], r["cut_left"], d)].add(r["sample"])
shared = {f"{k[0]}:{k[1]} {k[2]}": sorted(v) for k, v in rep.items() if len(v) > 1}

# ───────── ⑤ Sample_4 的 chr19 INV ─────────
inv = []
for reg, tag in [("chr19:46015900-46016400", "INV_breakpoints"), ("chr19:46016550-46016850", "offtarget_candidate")]:
    for smp in names:
        out = subprocess.run(["samtools", "view", "-c", "-F", "0x400", f"{BAM}/{smp}.markdup.bam", reg],
                             capture_output=True, text=True).stdout.strip()
        sup = subprocess.run(["samtools", "view", "-c", "-f", "0x800", f"{BAM}/{smp}.markdup.bam", reg],
                             capture_output=True, text=True).stdout.strip()
        inv.append(dict(region=reg, what=tag, sample=smp, group=role[smp],
                        reads=int(out or 0), supplementary=int(sup or 0)))
with open(f"{OUT}/tables/sample4_chr19_inv_context.tsv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(inv[0]), delimiter="\t"); w.writeheader(); w.writerows(inv)

# ───────── ③ 载体拷贝数按候选读段总数归一 ─────────
cand = {}
for smp in names:
    a = {}
    for ln in open(f"{S3}/{smp}.adapter_audit.tsv"):
        p = ln.split("\t")
        if len(p) == 2: a[p[0]] = int(p[1])
    cand[smp] = a.get("cand", 0)
med = sorted(cand.values())[len(cand) // 2]
cov = [r for r in csv.DictReader(open(f"{S3}/vector_coverage_n10.tsv"), delimiter="\t")]
norm = []
for r in cov:
    if r["region"] not in ("WHOLE_VECTOR", "AmpR", "Cas9", "EGFP", "mCherry", "ori"): continue
    s = r["sample"]; raw = r["est_copies_per_diploid"]
    try: raw = float(raw)
    except ValueError: continue
    norm.append(dict(sample=s, group=role[s], vector=r["vector"], region=r["region"],
                     est_copies_raw=round(raw, 3), candidate_reads=cand[s],
                     norm_factor=round(med / cand[s], 3) if cand[s] else "NA",
                     est_copies_normalized=round(raw * med / cand[s], 3) if cand[s] else "NA"))
with open(f"{OUT}/tables/vector_copies_normalized.tsv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(norm[0]), delimiter="\t"); w.writeheader(); w.writerows(norm)

# ───────── ④ Stage 1 定标表修正（原表不动）─────────
cal = list(csv.DictReader(open(f"{S1}/calibration_vs_client.tsv"), delimiter="\t"))
fixed, undet = [], 0
for r in cal:
    call = r["pipeline_call"]
    agrees = r["agrees"]
    if r["sample"] == "Sample_6":
        agrees, undet = "undetermined", 1
        note = "靶点窗口读段 19,620（其余样本 982–1,290）、d_norm=3.26（缺失区深度反而升高）⇒ amplicon 样污染，基因型不可按表面值判读"
    else:
        note = ""
    fixed.append(dict(sample=r["sample"], prior_from_client=r["prior_from_client"],
                      pipeline_call=call, agrees=agrees, note=note))
n_yes = sum(1 for r in fixed if r["agrees"] == "yes")
with open(f"{OUT}/tables/calibration_vs_client_corrected.tsv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(fixed[0]), delimiter="\t"); w.writeheader(); w.writerows(fixed)

summary = dict(
    item1_treated_only=dict(sites=len(sites), evidence_reads=len(rows),
                            same_indel_in_multiple_samples=shared,
                            verdict=("无任何 indel 在多个样本重复出现 ⇒ 符合各克隆独立本底，不支持真脱靶"
                                     if not shared else "存在跨样本重复的 indel ⇒ 需逐条判读")),
    item3_vector_norm=dict(candidate_reads_per_sample=cand, median_used=med),
    item4_calibration=dict(agrees_yes=n_yes, undetermined=undet,
                           corrected_statement=f"{n_yes}/10 与客户凝胶/Sanger 一致 + {undet} 个未定（Sample_6）"),
    note="原 calibration_vs_client.tsv 未改动；本目录为修正版，报告引用本目录。")
json.dump(summary, open(f"{OUT}/tables/verify_summary.json", "w"), indent=1, ensure_ascii=False)
print(json.dumps(summary, indent=1, ensure_ascii=False))
