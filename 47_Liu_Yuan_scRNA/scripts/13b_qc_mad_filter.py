#!/usr/bin/env python
"""P47 round1 step 13b：逐样本 MAD 自适应 QC + 并入 doublet + filter-impact audit。

上游：ckpt_12_raw_concat.h5ad (step 12) + 13a_doublet_percell.csv (step 13a)
产出：ckpt_13_qc.h5ad + 逐步过滤审计表 + machine-readable excluded table

判据：/scrna §3.2 —— 逐样本 MAD 自适应阈值，不用固定阈值(固定阈值已不被接受为默认,
且系统性丢低 RNA 含量细胞)。本批是单一细胞系 293T，细胞类型异质性低，但 MAD 仍逐样本算,
因为样本间测序深度差 2.98×，全局阈值会整体判掉低深度样本。
"""
import os, sys
for v in ("OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS","NUMEXPR_NUM_THREADS"):
    os.environ[v] = "16"
import numpy as np, pandas as pd, scanpy as sc

CKPT = "/mnt/ex_8T_SSD/projects_2026Q4/47_Liu_Yuan_scRNA/ckpt_1125"
OUT  = "/home/gao/projects_2026Q4/47_Liu_Yuan_scRNA/output_results/round1_descriptive_1125"
print("python :", sys.executable, flush=True)

src = f"{CKPT}/ckpt_12_raw_concat.h5ad"
assert os.path.exists(src), f"缺上游 checkpoint（由 step 12 产出）: {src}"
adata = sc.read_h5ad(src)
print(f"loaded {adata.shape}", flush=True)

dbl = pd.read_csv(f"{OUT}/13a_doublet_percell.csv").set_index("cell")
missing = adata.obs_names.difference(dbl.index)
assert len(missing) == 0, f"{len(missing)} 个细胞在 doublet 表里找不到（身份对齐失败）: {list(missing[:3])}"
# 按名字 reindex，绝不按位置赋值（按位置给输出补名字 = 静默错位，proj37 实测错 29% 行）
adata.obs["scDblFinder_score"] = dbl.reindex(adata.obs_names)["scDblFinder_score"].values
adata.obs["scDblFinder_class"] = dbl.reindex(adata.obs_names)["scDblFinder_class"].values
assert adata.obs["scDblFinder_class"].isna().sum() == 0

def mad_outlier(x, nmads=3.0, log=True):
    """isOutlier 等价：逐样本自适应。返回 (is_low, is_high, lo, hi)。"""
    v = np.log1p(x) if log else np.asarray(x, float)
    med = np.median(v); mad = np.median(np.abs(v - med))
    sigma = 1.4826 * mad
    lo, hi = med - nmads * sigma, med + nmads * sigma
    inv = (lambda z: np.expm1(z)) if log else (lambda z: z)
    return v < lo, v > hi, inv(lo), inv(hi)

audit, thr_rows = [], []
keep = pd.Series(True, index=adata.obs_names)
reason = pd.Series("", index=adata.obs_names)

for smp, idx in adata.obs.groupby("sample_id", observed=True).groups.items():
    o = adata.obs.loc[idx]
    n0 = len(o)
    lo_c, hi_c, loc_v, hic_v = mad_outlier(o.total_counts.values)
    lo_g, hi_g, log_v, hig_v = mad_outlier(o.n_genes_by_counts.values)
    # 线粒体比例：只切高端（低 mt 不是质量问题）；用 MAD 而非固定 5%
    _, hi_m, _, him_v = mad_outlier(o.pct_counts_mt.values, nmads=3.0, log=False)
    is_dbl = (o.scDblFinder_class == "doublet").values

    # ⚠ doublet 只记录、【不并入 bad】：293T 是单一细胞系，异型 doublet 原理上不存在，
    # scDblFinder 在纯细胞系里很可能把 S/G2M 高复杂度细胞判成 doublet（实测 14.8–22.3%）。
    # 硬删 19.6% 的细胞属于改变 data universe ⇒ 按 decision table §5 不自行决定。
    # 主口径 = 仅 MAD QC；另给 dbl_keep 列作敏感性臂，两条臂共用一个 checkpoint。
    flags = {"low_total_counts": lo_c, "high_total_counts": hi_c,
             "low_n_genes": lo_g, "high_n_genes": hi_g,
             "high_pct_mt": hi_m}
    bad = np.zeros(n0, bool)
    for name, m in flags.items():
        audit.append(dict(sample_id=smp, filter=name, n_flagged=int(m.sum()),
                          pct_flagged=round(100*float(m.mean()), 3), applied_in_primary=True))
        bad |= m
    audit.append(dict(sample_id=smp, filter="doublet_flag_not_applied",
                      n_flagged=int(is_dbl.sum()), pct_flagged=round(100*float(is_dbl.mean()), 3),
                      applied_in_primary=False))
    flags["doublet"] = is_dbl        # 仅用于下面记录首要原因
    # 记录首要原因（按优先级，便于人读；machine-readable 表里每条 flag 都单列）
    r = np.array([""]*n0, dtype=object)
    for name in ["high_pct_mt","low_n_genes","low_total_counts","high_n_genes","high_total_counts"]:
        r[(r == "") & flags[name]] = name
    keep.loc[o.index] = ~bad
    reason.loc[o.index] = r
    thr_rows.append(dict(sample_id=smp, n_before=n0, n_after=int((~bad).sum()),
                         n_removed=int(bad.sum()), pct_removed=round(100*float(bad.mean()),3),
                         total_counts_lo=round(float(loc_v),1), total_counts_hi=round(float(hic_v),1),
                         n_genes_lo=round(float(log_v),1), n_genes_hi=round(float(hig_v),1),
                         pct_mt_hi=round(float(him_v),3)))
    print(f"  {smp:12s} {n0:6d} -> {int((~bad).sum()):6d}  removed {int(bad.sum()):5d} ({100*bad.mean():5.2f}%)"
          f"  thr: counts[{loc_v:.0f},{hic_v:.0f}] genes[{log_v:.0f},{hig_v:.0f}] mt<={him_v:.2f}%", flush=True)

pd.DataFrame(audit).to_csv(f"{OUT}/13b_filter_impact_audit.csv", index=False)
pd.DataFrame(thr_rows).sort_values("sample_id").to_csv(f"{OUT}/13b_per_sample_thresholds.csv", index=False)
excl = adata.obs.loc[~keep.values, ["sample_id","pool","group","total_counts","n_genes_by_counts",
                                    "pct_counts_mt","scDblFinder_score","scDblFinder_class"]].copy()
excl["exclusion_reason"] = reason.loc[~keep.values].values
excl.to_csv(f"{OUT}/13b_excluded_cells.csv.gz", compression="gzip")
print(f"\nexcluded table: {excl.shape[0]} rows -> 13b_excluded_cells.csv.gz", flush=True)

adata.obs["qc_keep"] = keep.values
adata.obs["qc_exclusion_reason"] = reason.values
n_before = adata.n_obs
adata = adata[keep.values].copy()
print(f"\n主口径(仅 MAD QC) TOTAL {n_before} -> {adata.n_obs}  removed {n_before-adata.n_obs}"
      f" ({100*(n_before-adata.n_obs)/n_before:.2f}%)")
# 敏感性臂：在主口径之上再去 doublet。只加一列布尔，不另存 checkpoint。
adata.obs["dbl_keep"] = (adata.obs["scDblFinder_class"] != "doublet").values
arm = adata.obs.groupby("sample_id", observed=True).agg(
    n_primary=("dbl_keep", "size"), n_sens_no_doublet=("dbl_keep", "sum"))
arm["n_doublet_in_primary"] = arm.n_primary - arm.n_sens_no_doublet
arm["pct_doublet_in_primary"] = (100*arm.n_doublet_in_primary/arm.n_primary).round(2)
arm.sort_index().to_csv(f"{OUT}/13b_two_arms_cellcounts.csv")
print(f"敏感性臂(再去 doublet) -> {int(adata.obs.dbl_keep.sum())} cells"
      f" (再去掉 {int((~adata.obs.dbl_keep).sum())}, {100*(~adata.obs.dbl_keep).mean():.2f}%)")
print(arm.sort_index().to_string())

# 过滤后每样本仍需 >= 1000 细胞才有描述性价值（本批最小样本 12,027，门槛宽松）
post = adata.obs.sample_id.value_counts()
assert post.min() >= 1000, f"某样本过滤后细胞数过低: {post.to_dict()}"
# 基因侧：只去掉全 0 基因，并记录（HVG 只打标记，绝不子集化基因 —— 留到 step 14）
nz = np.asarray((adata.layers["counts"] > 0).sum(0)).ravel()
pd.DataFrame(dict(n_genes_total=[adata.n_vars], n_genes_zero_everywhere=[int((nz==0).sum())],
                  n_genes_kept=[int((nz>0).sum())])).to_csv(f"{OUT}/13b_gene_filter_audit.csv", index=False)
adata = adata[:, nz > 0].copy()
print(f"genes {len(nz)} -> {adata.n_vars} (去掉全 0 基因 {int((nz==0).sum())} 个)")

post.sort_index().to_frame("n_cells_after_qc").to_csv(f"{OUT}/13b_cells_after_qc.csv")
adata.write(f"{CKPT}/ckpt_13_qc.h5ad")
print(f"wrote {CKPT}/ckpt_13_qc.h5ad  shape={adata.shape}")
print("STEP13B_PASS")
