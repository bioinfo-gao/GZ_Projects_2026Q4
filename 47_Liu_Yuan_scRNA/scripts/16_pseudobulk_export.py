#!/usr/bin/env python
"""P47 round1 step 16：pseudobulk 聚合与导出（客户唯一能自行复分析的资产）。

上游：ckpt_14_clustered.h5ad   产出：sample 级与 sample×cluster 级 raw-count 矩阵 + coldata

判据：/scrna §3.5 —— 聚合的必须是 raw counts（整数），每个"细胞类型"单独建模，
      剔除 <10 细胞的格子并列出剔除了哪些。
⚠ 本项目不做跨条件建模（n≤2，见 decision table §1）⇒ 这里只导出资产，不建模、不出 p 值。
⚠ 本项目是单一细胞系 293T，"细胞类型"维度用 leiden 簇代替，且簇主要反映细胞周期
   ⇒ 簇级 pseudobulk 只作探索，不可解释为细胞类型特异效应。
"""
import os, sys, json
for v in ("OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS","NUMEXPR_NUM_THREADS"):
    os.environ[v] = "16"
import numpy as np, pandas as pd, scanpy as sc, scipy.sparse as sp

CKPT = "/mnt/ex_8T_SSD/projects_2026Q4/47_Liu_Yuan_scRNA/ckpt_1125"
OUT  = "/home/gao/projects_2026Q4/47_Liu_Yuan_scRNA/output_results/round1_descriptive_1125"
MIN_CELLS = 10
print("python :", sys.executable, flush=True)
src = f"{CKPT}/ckpt_14_clustered.h5ad"
assert os.path.exists(src), f"缺上游 checkpoint（由 step 14 产出）: {src}"
adata = sc.read_h5ad(src)
print(f"loaded {adata.shape}", flush=True)
C = adata.layers["counts"]
assert sp.issparse(C) and np.all(C.data == np.round(C.data)), "counts 不是整数稀疏矩阵"
# float32 求和超过 2^24 即失精度（样本总量 ~1e9，1125 首跑互证差 1–60 count）⇒ 全程 float64
C = C.astype(np.float64)

def agg(keys_df, label):
    """按 keys_df 的行分组求和 raw counts。返回 (counts_df[group x gene], coldata)。"""
    g = keys_df.astype(str).agg("|".join, axis=1)
    cats = pd.Index(sorted(g.unique()))
    ind = pd.Categorical(g, categories=cats).codes
    M = sp.csr_matrix((np.ones(len(ind)), (ind, np.arange(len(ind)))), shape=(len(cats), len(ind)))
    agg_counts = M @ C                                  # group x gene
    df = pd.DataFrame(np.asarray(agg_counts.todense()), index=cats, columns=adata.var_names).astype(np.int64)
    n = pd.Series(np.bincount(ind, minlength=len(cats)), index=cats, name="n_cells")
    return df, n

# ---------- (A) sample 级 ----------
samp_counts, samp_n = agg(adata.obs[["sample_id"]], "sample")
samp_counts.index.name = "sample_id"
# 独立互证：聚合后的总 count 必须等于逐细胞总和（拆包/聚合后必须重算，不能沿用）
tot_cells = np.asarray(C.sum(1)).ravel()
per_sample_direct = pd.Series(tot_cells, index=adata.obs_names).groupby(adata.obs.sample_id.values).sum()
chk = samp_counts.sum(1).astype(np.int64).sort_index()
assert (chk.values == per_sample_direct.sort_index().values.astype(np.int64)).all(), \
    "sample 级聚合总量与逐细胞直接求和不一致"
print("OK  sample 级聚合总量 == 逐细胞直接求和（独立互证通过）", flush=True)

coldata = adata.obs.groupby("sample_id", observed=True).agg(
    pool=("pool","first"), group=("group","first"), editor=("editor","first"), guide=("guide","first"),
    n_cells=("sample_id","size"), median_umi=("total_counts","median"),
    median_genes=("n_genes_by_counts","median"), mean_pct_mt=("pct_counts_mt","mean"),
    pct_S=("phase", lambda s: 100*float((s=="S").mean())),
    pct_G2M=("phase", lambda s: 100*float((s=="G2M").mean())),
    pct_doublet_flag=("scDblFinder_class", lambda s: 100*float((s=="doublet").mean()))).round(4)
coldata["total_pseudobulk_umi"] = chk
coldata.to_csv(f"{OUT}/16_pseudobulk_sample_coldata.csv")
samp_counts.T.to_csv(f"{OUT}/16_pseudobulk_sample_counts_genes_x_samples.csv.gz", compression="gzip")
print(f"wrote sample-level pseudobulk: {samp_counts.shape[1]} genes x {samp_counts.shape[0]} samples", flush=True)

# ---------- (B) sample x leiden 级，带 <10 细胞剔除清单 ----------
clus_counts, clus_n = agg(adata.obs[["sample_id","leiden"]], "sample_cluster")
kept = clus_n[clus_n >= MIN_CELLS].index
dropped = clus_n[clus_n < MIN_CELLS]
print(f"sample x cluster 格子: {len(clus_n)}  保留 {len(kept)}  剔除(<{MIN_CELLS} 细胞) {len(dropped)}", flush=True)
dropped.to_frame("n_cells").to_csv(f"{OUT}/16_pseudobulk_dropped_cells_lt10.csv")
clus_counts.loc[kept].T.to_csv(f"{OUT}/16_pseudobulk_sample_x_cluster_counts.csv.gz", compression="gzip")
cl_col = pd.DataFrame({"sample_cluster": kept})
cl_col[["sample_id","leiden"]] = cl_col.sample_cluster.str.split("|", expand=True)
cl_col = cl_col.merge(coldata[["pool","group","editor","guide"]], left_on="sample_id", right_index=True)
cl_col["n_cells"] = clus_n.loc[kept].values
cl_col.to_csv(f"{OUT}/16_pseudobulk_sample_x_cluster_coldata.csv", index=False)

json.dump({"min_cells_per_cell_in_pseudobulk": MIN_CELLS,
           "n_samples": int(samp_counts.shape[0]), "n_genes": int(samp_counts.shape[1]),
           "n_sample_cluster_cells_total": int(len(clus_n)),
           "n_sample_cluster_cells_kept": int(len(kept)),
           "n_sample_cluster_cells_dropped": int(len(dropped)),
           "no_modeling_performed": True,
           "reason_no_pvalues": "n<=2 per condition; pool fully confounded with group; see decision table"},
          open(f"{OUT}/16_pseudobulk_facts.json","w"), indent=2)
print(coldata.to_string())
print("STEP16_PASS")
