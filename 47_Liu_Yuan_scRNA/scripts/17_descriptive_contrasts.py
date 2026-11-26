#!/usr/bin/env python
"""P47 round1 step 17：描述性效应量（**不出 p 值 / 不出 FDR**）。

上游：16_pseudobulk_sample_counts_genes_x_samples.csv.gz + coldata
产出：各对比的 log2FC 表 + "pool 效应 vs pool 内效应"的量级对照

为什么不出 p 值（三条独立理由，任一条都足够，见 decision table §1）：
  1. n<=2/条件，低于 /scrna §1.1 的 n>=3 跨条件推断门槛；
  2. edited/control 与 pool 完全共线，主方向不可识别；
  3. Flex 同 pool 内 GEM 重叠 11.4-17.4%，样本并非完全独立观测。
⇒ 本步只给方向与幅度。结论里不得出现"显著"、也不得出现"无变化"。

本步的核心产物不是基因清单，而是最后那张量级对照表：
  如果【跨 pool 的混杂对比】检出的变化远多于【pool 内可识别对比】，
  那就直接用数据证明了主对比里看到的差异主要是 pool/批次效应。
"""
import os, sys, json
import numpy as np, pandas as pd

OUT = "/home/gao/projects_2026Q4/47_Liu_Yuan_scRNA/output_results/round1_descriptive_1125"
print("python :", sys.executable, flush=True)
cnt = pd.read_csv(f"{OUT}/16_pseudobulk_sample_counts_genes_x_samples.csv.gz", index_col=0)
col = pd.read_csv(f"{OUT}/16_pseudobulk_sample_coldata.csv", index_col=0)
cnt = cnt[col.index]                                    # 列序对齐（按名字，不按位置）
assert list(cnt.columns) == list(col.index)
print(f"pseudobulk: {cnt.shape[0]} genes x {cnt.shape[1]} samples", flush=True)

# CPM（不是 DESeq2 size factor：本步刻意不建模，避免顺手产出 p 值）
cpm = cnt / cnt.sum(0) * 1e6
EXPR_MIN = 1.0                                          # CPM>=1 的基因才参与（低表达的 logFC 噪声大）
expressed = (cpm >= EXPR_MIN).sum(1) >= 4
print(f"表达过滤: CPM>={EXPR_MIN} 于 >=4 个样本 -> {int(expressed.sum())}/{len(expressed)} 基因参与", flush=True)
lg = np.log2(cpm.loc[expressed] + 1)

CONTRASTS = {
  "pool1_ABE_vs_CBE":       (["N1-ABE","N2-ABE"],      ["N1-CBE","N2-CBE"],          "pool1", "identifiable_within_pool"),
  "pool1_guideN1_vs_N2":    (["N1-ABE","N1-CBE"],      ["N2-ABE","N2-CBE"],          "pool1", "identifiable_within_pool"),
  "pool2_editoronly_vs_gRNAonly": (["ABE-Control","CBE-Control"], ["N1-Control","N2-Control"], "pool2", "identifiable_within_pool"),
  "pool2_ABEctrl_vs_CBEctrl":(["ABE-Control"],         ["CBE-Control"],              "pool2", "identifiable_within_pool_n1"),
  "CONFOUNDED_edited_vs_control": (["N1-ABE","N2-ABE","N1-CBE","N2-CBE"],
                                   ["N1-Control","N2-Control","CBE-Control","ABE-Control"],
                                   "cross-pool", "CONFOUNDED_with_pool_not_interpretable"),
}
rows, summary = {}, []
for name, (a, b, scope, status) in CONTRASTS.items():
    fc = lg[a].mean(1) - lg[b].mean(1)
    df = pd.DataFrame({"log2FC": fc.round(4),
                       "mean_log2CPM_groupA": lg[a].mean(1).round(4),
                       "mean_log2CPM_groupB": lg[b].mean(1).round(4),
                       "mean_log2CPM_all": lg.mean(1).round(4)})
    df["abs_log2FC"] = df.log2FC.abs()
    df = df.sort_values("abs_log2FC", ascending=False)
    df.index.name = "gene"
    df.to_csv(f"{OUT}/17_contrast_{name}.csv.gz", compression="gzip")
    rows[name] = df
    summary.append(dict(contrast=name, scope=scope, status=status,
                        n_groupA=len(a), n_groupB=len(b),
                        n_genes_tested=int(df.shape[0]),
                        n_abs_log2FC_gt_0_5=int((df.abs_log2FC > 0.5).sum()),
                        n_abs_log2FC_gt_1=int((df.abs_log2FC > 1).sum()),
                        n_abs_log2FC_gt_2=int((df.abs_log2FC > 2).sum()),
                        median_abs_log2FC=round(float(df.abs_log2FC.median()), 4),
                        p95_abs_log2FC=round(float(df.abs_log2FC.quantile(.95)), 4)))
    print(f"  {name:34s} |log2FC|>1: {int((df.abs_log2FC>1).sum()):5d}  "
          f"median|log2FC|={df.abs_log2FC.median():.4f}  p95={df.abs_log2FC.quantile(.95):.4f}", flush=True)

sm = pd.DataFrame(summary)
sm.to_csv(f"{OUT}/17_contrast_summary.csv", index=False)

# ---------- 技术重复级本底：同 pool 内、同处理类别的两个样本互比 ----------
# 这是"这条流程自己的噪声本底"：N1-ABE vs N2-ABE 只差 gRNA，N1-Control vs N2-Control 同理。
BASELINE = {"within_pool1_ABE_N1_vs_N2": ("N1-ABE","N2-ABE"),
            "within_pool1_CBE_N1_vs_N2": ("N1-CBE","N2-CBE"),
            "within_pool2_gRNActrl_N1_vs_N2": ("N1-Control","N2-Control")}
base_rows = []
for name, (x, y) in BASELINE.items():
    fc = (lg[x] - lg[y]).abs()
    base_rows.append(dict(pair=name, n_abs_log2FC_gt_1=int((fc > 1).sum()),
                          median_abs_log2FC=round(float(fc.median()), 4),
                          p95_abs_log2FC=round(float(fc.quantile(.95)), 4)))
    print(f"  [本底] {name:32s} |log2FC|>1: {int((fc>1).sum()):5d}  median={fc.median():.4f}", flush=True)
bl = pd.DataFrame(base_rows); bl.to_csv(f"{OUT}/17_single_sample_pair_baseline.csv", index=False)

# ---------- 核心对照：混杂对比 vs pool 内对比 vs 单样本本底 ----------
conf = sm.loc[sm.contrast == "CONFOUNDED_edited_vs_control", "n_abs_log2FC_gt_1"].iloc[0]
within = sm.loc[sm.status == "identifiable_within_pool", "n_abs_log2FC_gt_1"]
base_med = int(bl.n_abs_log2FC_gt_1.median())
verdict = pd.DataFrame([
    dict(comparison="CONFOUNDED cross-pool edited_vs_control", n_genes_abs_log2FC_gt_1=int(conf)),
    dict(comparison="identifiable within-pool contrasts (median)", n_genes_abs_log2FC_gt_1=int(within.median())),
    dict(comparison="single-sample-pair technical-ish baseline (median)", n_genes_abs_log2FC_gt_1=base_med),
])
verdict["ratio_vs_within_pool_median"] = (verdict.n_genes_abs_log2FC_gt_1 /
                                          max(int(within.median()), 1)).round(2)
verdict.to_csv(f"{OUT}/17_pool_effect_vs_within_pool_effect.csv", index=False)
print("\n### 量级对照（本步核心产物）：\n" + verdict.to_string(index=False), flush=True)

json.dump({"no_pvalues_reported": True,
           "reasons": ["n<=2 per condition (scrna SOP n>=3 gate)",
                       "group fully collinear with GEX pool",
                       "Flex within-pool GEM overlap 11.4-17.4% breaks independence"],
           "expr_filter": f"CPM>={EXPR_MIN} in >=4 samples",
           "n_genes_tested": int(lg.shape[0]),
           "confounded_contrast_n_gt1": int(conf),
           "within_pool_median_n_gt1": int(within.median()),
           "single_pair_baseline_median_n_gt1": base_med},
          open(f"{OUT}/17_contrast_facts.json","w"), indent=2)
print("STEP17_PASS")
