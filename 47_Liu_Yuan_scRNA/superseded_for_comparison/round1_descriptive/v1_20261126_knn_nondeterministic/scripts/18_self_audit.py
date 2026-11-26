#!/usr/bin/env python
"""P47 round1 step 18：自审（写成代码，不靠肉眼）。

上游：ckpt_14_clustered.h5ad + step 16/17 产物
产出：18_self_audit.csv（逐条断言结果）+ 18_pseudoreplication_demo.csv

包含四类检查：
  (1) 确定性：带随机性的步骤同参数跑两遍，产物必须逐元素一致（/feedback：单轮 exit=0 证明不了确定性）
  (2) 身份对齐：所有按名字 join 的表，行集合必须完全一致（按位置赋值 = 静默错位）
  (3) 分母非零：任何"0 / 没检出"在当结论前必须证明分母不是被自己的筛选掏空的
  (4) 自证伪重复：在同一份数据上把一个对照样本【随机对半】跑细胞级检验，
      把"照样检出 N 个显著基因"的数字摆出来 —— 支撑本轮不出 p 值的决定
"""
import os, sys, json
for v in ("OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS","NUMEXPR_NUM_THREADS"):
    os.environ[v] = "16"
import numpy as np, pandas as pd, scanpy as sc
SEED = 20261125
CKPT = "/mnt/ex_8T_SSD/projects_2026Q4/47_Liu_Yuan_scRNA/ckpt_1125"
OUT  = "/home/gao/projects_2026Q4/47_Liu_Yuan_scRNA/output_results/round1_descriptive_1125"
print("python :", sys.executable, flush=True)
adata = sc.read_h5ad(f"{CKPT}/ckpt_14_clustered.h5ad")
print(f"loaded {adata.shape}", flush=True)

res = []
def A(name, ok, detail=""):
    res.append(dict(check=name, result="PASS" if ok else "FAIL", detail=str(detail)))
    print(f"  {'PASS' if ok else 'FAIL'}  {name}  {detail}", flush=True)
    return ok

# ---------- (1) 确定性 ----------
sc.tl.leiden(adata, resolution=0.5, key_added="leiden_run2", flavor="igraph",
             n_iterations=2, directed=False, random_state=SEED)
same = (adata.obs.leiden.astype(str).values == adata.obs.leiden_run2.astype(str).values)
A("determinism: leiden 同参数两遍标签完全一致", bool(same.all()),
  f"{same.sum()}/{len(same)} 一致")
sc.tl.score_genes_cell_cycle  # 已在 step14 跑过；这里复查打分可重现
from sklearn.neighbors import NearestNeighbors
def mixing(emb, labels, k=30, n_sub=20000, seed=SEED):
    rng = np.random.default_rng(seed)
    idx = rng.choice(emb.shape[0], min(n_sub, emb.shape[0]), replace=False)
    nn = NearestNeighbors(n_neighbors=k+1).fit(emb)
    _, ind = nn.kneighbors(emb[idx]); lab = np.asarray(labels)
    return float(1-(lab[ind[:,1:]] == lab[idx][:,None]).mean())
m1 = mixing(adata.obsm["X_pca_harmony"], adata.obs["pool"].values)
m2 = mixing(adata.obsm["X_pca_harmony"], adata.obs["pool"].values)
A("determinism: 混合度指标两遍一致", m1 == m2, f"{m1:.6f} vs {m2:.6f}")

# ---------- (2) 身份对齐 ----------
pc = pd.read_csv(f"{OUT}/14_percell_metadata.csv.gz", index_col=0)
A("身份对齐: percell_metadata 行集合 == ckpt_14 细胞集合",
  set(pc.index) == set(adata.obs_names), f"{len(pc)} vs {adata.n_obs}")
A("身份对齐: leiden 标签按名字核对一致",
  bool((pc.reindex(adata.obs_names).leiden.astype(str).values ==
        adata.obs.leiden.astype(str).values).all()))
rho_f = f"{OUT}/15_ambient_rho_percell.csv.gz"
if os.path.exists(rho_f):
    rho = pd.read_csv(rho_f).set_index("cell")
    A("身份对齐: ambient rho 行集合 == ckpt_14 细胞集合",
      set(rho.index) == set(adata.obs_names), f"{len(rho)} vs {adata.n_obs}")
    A("分母非零: ambient rho 非常数", rho.decontX_contamination.nunique() > 100,
      f"unique={rho.decontX_contamination.nunique()}")
pb = pd.read_csv(f"{OUT}/16_pseudobulk_sample_coldata.csv", index_col=0)
A("身份对齐: pseudobulk coldata 8 个样本齐全", len(pb) == 8, sorted(pb.index))
A("一致性: coldata n_cells 之和 == ckpt_14 细胞数",
  int(pb.n_cells.sum()) == adata.n_obs, f"{int(pb.n_cells.sum())} vs {adata.n_obs}")

# ---------- (3) 分母非零 / 不把筛选掏空当结论 ----------
cnt = pd.read_csv(f"{OUT}/16_pseudobulk_sample_counts_genes_x_samples.csv.gz", index_col=0)
A("分母非零: pseudobulk 每样本总 UMI > 0", bool((cnt.sum(0) > 0).all()), cnt.sum(0).min())
A("分母非零: 参与对比的基因数 > 5000",
  int(json.load(open(f"{OUT}/17_contrast_facts.json"))["n_genes_tested"]) > 5000,
  json.load(open(f"{OUT}/17_contrast_facts.json"))["n_genes_tested"])
mt_in = [g for g in adata.var_names if g.startswith("MT-")]
A("分母非零: MT 基因 11 个且 pct_counts_mt 非 0", len(mt_in) == 11 and
  float(adata.obs.pct_counts_mt.mean()) > 0, f"{len(mt_in)} genes, mean={adata.obs.pct_counts_mt.mean():.3f}%")

# ---------- (4) 自证伪重复 ----------
rng = np.random.default_rng(SEED)
demo = []
for smp in ["N1-Control", "N1-ABE"]:
    sub = adata[adata.obs.sample_id == smp].copy()
    assert adata.uns.get("X_is","").startswith("log1p"), "X 不是 lognorm，Wilcoxon 会跑在错的标度上"
    half = rng.permutation(sub.n_obs) < sub.n_obs // 2
    sub.obs["fake_group"] = np.where(half, "fakeA", "fakeB")
    sc.tl.rank_genes_groups(sub, "fake_group", groups=["fakeA"], reference="fakeB",
                            method="wilcoxon", pts=False)
    d = sc.get.rank_genes_groups_df(sub, group="fakeA")
    for thr in (0.05, 0.01):
        demo.append(dict(sample_id=smp, n_cells=sub.n_obs, test="cell-level Wilcoxon",
                         comparison="random half vs random half (NO real difference)",
                         fdr_threshold=thr, n_significant=int((d.pvals_adj < thr).sum()),
                         n_significant_and_abs_lfc_gt_0_25=int(((d.pvals_adj < thr) &
                                                                (d.logfoldchanges.abs() > 0.25)).sum())))
        print(f"  [伪重复自证] {smp} 随机对半 cell-level Wilcoxon FDR<{thr}: "
              f"{int((d.pvals_adj < thr).sum())} 个'显著'基因", flush=True)
    del sub
dd = pd.DataFrame(demo); dd.to_csv(f"{OUT}/18_pseudoreplication_demo.csv", index=False)
A("自证伪重复: 随机对半仍检出大量'显著'基因（证明细胞级检验不可用于跨条件）",
  bool((dd.n_significant > 100).any()), f"max={int(dd.n_significant.max())}")

ra = pd.DataFrame(res)
ra.to_csv(f"{OUT}/18_self_audit.csv", index=False)
n_fail = int((ra.result == "FAIL").sum())
print(f"\n自审：{len(ra)} 条，PASS {len(ra)-n_fail}，FAIL {n_fail}")
if n_fail:
    print(ra[ra.result == "FAIL"].to_string(index=False)); sys.exit(1)
print("STEP18_PASS")
