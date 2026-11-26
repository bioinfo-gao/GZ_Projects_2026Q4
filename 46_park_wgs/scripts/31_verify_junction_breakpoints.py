"""P46 交付前复核 ② —— 130 个非共有整合接头断点是否"在多条读段间固定"。

为什么要查：真整合接头的基因组断点与载体断点都是**同一个碱基**（病毒整合是一次事件，
此后克隆扩增只是复制它）⇒ 支持它的多条读段应指向同一个载体坐标。
随机本底（比对歧义、嵌合读段）则载体侧坐标分散。
判据（写死）：该接头支持读段的载体坐标跨度 `vec_pos_max - vec_pos_min`
  ≤ 50 bp ⇒ fixed（可信）；> 150 bp ⇒ scattered（疑为本底）；之间 ⇒ ambiguous。
另按"携带样本数"分层：10 样本共有 = 亲本背景（阳性对照所在层），
仅 1 个样本 = 最可能是本底或真正的克隆特异事件 —— 两者靠固定性区分。
"""
import csv, json, collections, os

PROJ = "/home/gao/projects_2026Q4/46_park_wgs"
S3 = f"{PROJ}/analysis/stage3_vector_v2_20261124/tables"
OUT = f"{PROJ}/analysis/stage5_verify_20261126/tables"
os.makedirs(OUT, exist_ok=True)
FIXED, SCATTER = 50, 150

J = list(csv.DictReader(open(f"{S3}/integration_junctions_n10.tsv"), delimiter="\t"))
share = {(r["chrom"], r["breakpoint"], r["side"]): int(r["n_samples"])
         for r in csv.DictReader(open(f"{S3}/integration_sharing_n10.tsv"), delimiter="\t")}
PC = {("chr11", "61804435", "L"), ("chr11", "61804455", "R"),
      ("chr4", "100591812", "L"), ("chr4", "100591830", "R")}

rows = []
for r in J:
    key = (r["chrom"], r["breakpoint"], r["side"])
    spread = int(r["vec_pos_max"]) - int(r["vec_pos_min"])
    cls = "fixed" if spread <= FIXED else ("scattered" if spread > SCATTER else "ambiguous")
    rows.append(dict(sample=r["sample"], group=r["group"], chrom=r["chrom"], breakpoint=r["breakpoint"],
                     side=r["side"], reads=int(r["reads"]), vector=r["vector"],
                     vec_span_bp=spread, fixedness=cls,
                     n_samples_sharing=share.get(key, 1),
                     is_positive_control=int(key in PC), gene=r["gene"]))
with open(f"{OUT}/junction_breakpoint_fixedness.tsv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0]), delimiter="\t"); w.writeheader(); w.writerows(rows)

def tab(sel, label):
    sub = [r for r in rows if sel(r)]
    c = collections.Counter(r["fixedness"] for r in sub)
    med = sorted(r["vec_span_bp"] for r in sub)
    return dict(label=label, n=len(sub), fixed=c["fixed"], ambiguous=c["ambiguous"], scattered=c["scattered"],
                median_span_bp=med[len(med) // 2] if med else None,
                median_reads=sorted(r["reads"] for r in sub)[len(sub) // 2] if sub else None)

layers = [tab(lambda r: r["is_positive_control"] == 1, "阳性对照（P31 两个共有接头）"),
          tab(lambda r: r["n_samples_sharing"] == 10, "10 样本共有（亲本背景）"),
          tab(lambda r: 2 <= r["n_samples_sharing"] <= 9, "部分样本共有（2–9）"),
          tab(lambda r: r["n_samples_sharing"] == 1, "仅 1 个样本")]

# 仅处理组且 fixed 的位点 —— 唯一可能被读成"编辑相关新整合"的一类
cand = [r for r in rows if r["n_samples_sharing"] == 1 and r["fixedness"] == "fixed" and r["group"] == "treated"]
summary = dict(total_junction_records=len(rows), layers=layers,
               treated_only_fixed=[{k: r[k] for k in ("sample", "chrom", "breakpoint", "side", "reads", "vec_span_bp", "gene")} for r in cand],
               rule=dict(fixed_le_bp=FIXED, scattered_gt_bp=SCATTER),
               reading="阳性对照层若以 fixed 为主、而『仅 1 个样本』层以 scattered 为主 ⇒ 判据可用，"
                       "且单样本接头多为本底；treated_only_fixed 是唯一需要逐条判读的清单。")
json.dump(summary, open(f"{OUT}/junction_fixedness_summary.json", "w"), indent=1, ensure_ascii=False)
print(json.dumps(summary, indent=1, ensure_ascii=False))
