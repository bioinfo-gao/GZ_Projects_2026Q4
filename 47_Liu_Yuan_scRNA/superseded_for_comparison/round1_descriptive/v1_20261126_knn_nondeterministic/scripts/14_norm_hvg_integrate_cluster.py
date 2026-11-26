#!/usr/bin/env python
"""P47 round1 step 14：归一化 → HVG(只打标记) → PCA → Harmony(按 pool) → 聚类 → 细胞周期。

上游：ckpt_13_qc.h5ad   产出：ckpt_14_clustered.h5ad + 整合/聚类/周期审计表

设计说明（为什么这么选，不是默认就好）：
  * 批次变量 = pool。两个 pool 是不同 physical GEX library，且组别与 pool 完全共线
    ⇒ 整合只为"看结构"，【绝不能】用整合后的嵌入去论证组间差异（会把处理效应一起抹掉）。
    因此同时保存未整合 PCA 与整合后嵌入，并对照检查。
  * 293T 是单一细胞系 ⇒ 不做细胞类型注释。预期主轴是细胞周期，故显式打分验证。
  * HVG 只打标记不子集化（/scrna 规则 5）；counts 已在 step 12 存入 layers['counts']。
  * 不用 sc.pp.scale：本矩阵稠密度 44.4%(nnz≈1.09e9)，全矩阵稠密化会 OOM。
    PCA 只在 HVG 子空间内做 zero-center，由 scanpy 的 mask_var 内部处理。
"""
import os, sys, json
for v in ("OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS","NUMEXPR_NUM_THREADS"):
    os.environ[v] = "16"
import numpy as np, pandas as pd, scanpy as sc
SEED = 20261125
np.random.seed(SEED)
sc.settings.n_jobs = 16

CKPT = "/mnt/ex_8T_SSD/projects_2026Q4/47_Liu_Yuan_scRNA/ckpt_1125"
OUT  = "/home/gao/projects_2026Q4/47_Liu_Yuan_scRNA/output_results/round1_descriptive_1125"
print("python :", sys.executable, flush=True)
src = f"{CKPT}/ckpt_13_qc.h5ad"
assert os.path.exists(src), f"缺上游 checkpoint（由 step 13b 产出）: {src}"
adata = sc.read_h5ad(src)
print(f"loaded {adata.shape}", flush=True)
assert "counts" in adata.layers, "layers['counts'] 丢了 —— 下游一切 DE/pseudobulk 失效"

# ---------- 归一化 ----------
adata.X = adata.layers["counts"].copy()
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
# 不另存 lognorm 层：X 此刻【就是】log-normalized，再存一份等于白占 8.7 GB
# (本矩阵 nnz≈1.09e9，每层约 8.7 GB)。下游要 lognorm 直接用 X，要 counts 用 layers['counts']。
adata.uns["X_is"] = "log1p(normalize_total(counts, 1e4))"
print("normalized (shifted-log, target_sum=1e4); X is lognorm, layers['counts'] is raw", flush=True)

# ---------- HVG：只打标记 ----------
sc.pp.highly_variable_genes(adata, n_top_genes=2000, batch_key="pool")
n_hvg = int(adata.var.highly_variable.sum())
print(f"HVG flagged: {n_hvg} (未子集化，var 维度仍为 {adata.n_vars})", flush=True)
assert n_hvg >= 1500, f"HVG 数异常: {n_hvg}"
assert adata.n_vars == adata.layers["counts"].shape[1], "基因维度被意外子集化"

# ---------- 细胞周期打分（Tirosh et al. 2016 标准列表）----------
S_GENES = """MCM5 PCNA TYMS FEN1 MCM2 MCM4 RRM1 UNG GINS2 MCM6 CDCA7 DTL PRIM1 UHRF1 MLF1IP
HELLS RFC2 RPA2 NASP RAD51AP1 GMNN WDR76 SLBP CCNE2 UBR7 POLD3 MSH2 ATAD2 RAD51 RRM2 CDC45 CDC6
EXO1 TIPIN DSCC1 BLM CASP8AP2 USP1 CLSPN POLA1 CHAF1B BRIP1 E2F8""".split()
G2M_GENES = """HMGB2 CDK1 NUSAP1 UBE2C BIRC5 TPX2 TOP2A NDC80 CKS2 NUF2 CKS1B MKI67 TMPO CENPF TACC3
FAM64A SMC4 CCNB2 CKAP2L CKAP2 AURKB BUB1 KIF11 ANP32E TUBB4B GTSE1 KIF20B HJURP CDCA3 HN1 CDC20
TTK CDC25C KIF2C RANGAP1 NCAPD2 DLGAP5 CDCA2 CDCA8 ECT2 KIF23 HMMR AURKA PSRC1 ANLN LBR CKAP5
CENPE CTCF NEK2 G2E3 GAS2L3 CBX5 CENPA""".split()
s_in  = [g for g in S_GENES   if g in adata.var_names]
g2_in = [g for g in G2M_GENES if g in adata.var_names]
print(f"cell cycle genes found: S {len(s_in)}/{len(S_GENES)}  G2M {len(g2_in)}/{len(G2M_GENES)}", flush=True)
# "0 值当结论前先证明分母非零"：基因找不到时打分恒为 0，长得像"没有周期结构"
assert len(s_in) >= 25 and len(g2_in) >= 30, f"周期基因命中过少(S={len(s_in)},G2M={len(g2_in)})，打分不可信"
sc.tl.score_genes_cell_cycle(adata, s_genes=s_in, g2m_genes=g2_in, random_state=SEED)
phase = adata.obs.groupby(["sample_id","phase"], observed=True).size().unstack(fill_value=0)
phase_pct = (100*phase.div(phase.sum(1), axis=0)).round(2)
phase.to_csv(f"{OUT}/14_cellcycle_phase_counts.csv"); phase_pct.to_csv(f"{OUT}/14_cellcycle_phase_pct.csv")
print("\n### 细胞周期相(%)：\n" + phase_pct.to_string(), flush=True)
assert (phase_pct[[c for c in ("S","G2M") if c in phase_pct]].sum(1) > 10).all(), \
    "某样本 S+G2M < 10% —— 293T 应明显 cycling，先怀疑读入/归一化出错"

# ---------- PCA（未整合）----------
sc.pp.pca(adata, n_comps=50, mask_var="highly_variable", svd_solver="arpack", random_state=SEED)
adata.obsm["X_pca_unintegrated"] = adata.obsm["X_pca"].copy()
print("PCA done (50 comps, HVG subspace)", flush=True)

# ---------- Harmony（批次 = pool）----------
# 直接调 harmonypy：2.0.0 的 Z_corr 已是 cells×PCs，scanpy 1.11.5 包装器再 .T 一次会形状报错（1125 实测）
import harmonypy
_ho = harmonypy.run_harmony(adata.obsm["X_pca"], adata.obs, ["pool"], random_state=SEED)
_Z = np.asarray(_ho.Z_corr)
if _Z.shape == (adata.obsm["X_pca"].shape[1], adata.n_obs):   # 旧版 harmonypy 方向
    _Z = _Z.T
assert _Z.shape == adata.obsm["X_pca"].shape, f"Harmony 输出形状 {_Z.shape} != PCA {adata.obsm['X_pca'].shape}"
adata.obsm["X_pca_harmony"] = _Z
print(f"Harmony done  shape={_Z.shape}", flush=True)

# ---------- 量化对照：整合前后，批次混合 vs 生物保留 ----------
def mixing(emb, labels, k=30, n_sub=20000, seed=SEED):
    """每个细胞的 k 近邻里"异批次比例"的均值；越接近理论值=混合越好。"""
    from sklearn.neighbors import NearestNeighbors
    rng = np.random.default_rng(seed)
    idx = rng.choice(emb.shape[0], min(n_sub, emb.shape[0]), replace=False)
    nn = NearestNeighbors(n_neighbors=k+1).fit(emb)
    _, ind = nn.kneighbors(emb[idx])
    lab = np.asarray(labels)
    same = (lab[ind[:,1:]] == lab[idx][:,None]).mean()
    return float(1-same)
pool_lab = adata.obs["pool"].values
prop = pd.Series(pool_lab).value_counts(normalize=True)
expect_diff = float(1 - (prop**2).sum())      # 完全混合时的理论异批次比例
rows = []
for name, key in [("unintegrated","X_pca_unintegrated"), ("harmony","X_pca_harmony")]:
    mx_pool = mixing(adata.obsm[key], pool_lab)
    mx_phase = mixing(adata.obsm[key], adata.obs["phase"].values)
    rows.append(dict(embedding=name, cross_pool_knn_frac=round(mx_pool,4),
                     cross_pool_expected_if_mixed=round(expect_diff,4),
                     batch_mixing_ratio=round(mx_pool/expect_diff,4),
                     cross_phase_knn_frac=round(mx_phase,4),
                     note="cross_phase 越低=细胞周期结构保留越好"))
    print(f"  {name:14s} cross-pool kNN={mx_pool:.4f} (expect {expect_diff:.4f}, ratio {mx_pool/expect_diff:.3f})"
          f"  cross-phase kNN={mx_phase:.4f}", flush=True)
pd.DataFrame(rows).to_csv(f"{OUT}/14_integration_metrics.csv", index=False)

# ---------- 聚类 + UMAP（仅展示）----------
sc.pp.neighbors(adata, n_neighbors=15, use_rep="X_pca_harmony", random_state=SEED)
sc.tl.leiden(adata, resolution=0.5, key_added="leiden", flavor="igraph",
             n_iterations=2, directed=False, random_state=SEED)
sc.tl.umap(adata, random_state=SEED)
nclus = adata.obs.leiden.nunique()
print(f"\nleiden clusters: {nclus}", flush=True)

ct = adata.obs.groupby(["leiden","phase"], observed=True).size().unstack(fill_value=0)
ct_pct = (100*ct.div(ct.sum(1), axis=0)).round(2)
ct_pct.to_csv(f"{OUT}/14_cluster_by_phase_pct.csv")
cs = adata.obs.groupby(["sample_id","leiden"], observed=True).size().unstack(fill_value=0)
cs.to_csv(f"{OUT}/14_cluster_by_sample_counts.csv")
(100*cs.div(cs.sum(1), axis=0)).round(2).to_csv(f"{OUT}/14_cluster_by_sample_pct.csv")
cp = adata.obs.groupby(["leiden","pool"], observed=True).size().unstack(fill_value=0)
cp_pct = (100*cp.div(cp.sum(1), axis=0)).round(2)
cp_pct.to_csv(f"{OUT}/14_cluster_by_pool_pct.csv")
print("\n### 每簇的细胞周期相组成(%)：\n" + ct_pct.to_string(), flush=True)
print("\n### 每簇的 pool 组成(%)：\n" + cp_pct.to_string(), flush=True)

json.dump({"n_hvg": n_hvg, "n_clusters": int(nclus), "seed": SEED,
           "s_genes_found": len(s_in), "g2m_genes_found": len(g2_in),
           "cross_pool_expected_if_mixed": expect_diff,
           "integration_metrics": rows,
           "caveat": "pool 与 group 完全共线；整合后嵌入不得用于论证组间差异"},
          open(f"{OUT}/14_cluster_facts.json","w"), indent=2)
# 逐细胞元数据导出：R 侧(step 15 ambient)唯一的桥。桥是文件不是 rpy2(/scrna §3.5)。
obs_cols = ["sample_id","pool","group","editor","guide","leiden","phase","S_score","G2M_score",
            "total_counts","n_genes_by_counts","pct_counts_mt","scDblFinder_score",
            "scDblFinder_class","dbl_keep"]
adata.obs[[c for c in obs_cols if c in adata.obs]].to_csv(f"{OUT}/14_percell_metadata.csv.gz",
                                                          compression="gzip")
print(f"wrote {OUT}/14_percell_metadata.csv.gz  rows={adata.n_obs}", flush=True)

adata.write(f"{CKPT}/ckpt_14_clustered.h5ad")
print(f"\nwrote {CKPT}/ckpt_14_clustered.h5ad  shape={adata.shape}")
print("STEP14_PASS")
