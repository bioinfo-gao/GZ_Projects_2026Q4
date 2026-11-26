#!/usr/bin/env python
"""P47 round1 step 19：图件（PDF only）。

上游：ckpt_14_clustered.h5ad + step 15/16/17 产物   产出：figures_1125/*.pdf
house style 与配色一律走 /scientific-figure-making 的 assets/pubfig.py，不自造配色。
⚠ 同图内语义不同的变量绝不共用配色：pool / phase / sample / leiden 各用不同 family。
⚠ 散点必须断言实际可见点数 == 期望 N（NaN 坐标或 shape 超标记数会静默丢点）。
"""
import os, sys
for v in ("OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS","NUMEXPR_NUM_THREADS"):
    os.environ[v] = "16"
sys.path.insert(0, "/home/gao/.claude/skills/scientific-figure-making/assets")
import numpy as np, pandas as pd, scanpy as sc, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pubfig
from pubfig import (PALETTE, PALETTE_FAMILIES, apply_publication_style, create_subplots,
                    finalize_figure, assert_no_palette_collision, assert_visible_on_white, contrast_vs_white)

CKPT = "/mnt/ex_8T_SSD/projects_2026Q4/47_Liu_Yuan_scRNA/ckpt_1125"
OUT  = "/home/gao/projects_2026Q4/47_Liu_Yuan_scRNA/output_results/round1_descriptive_1125"
FIG  = f"{OUT}/figures_1125"
os.makedirs(FIG, exist_ok=True)
apply_publication_style()
print("python :", sys.executable, flush=True)
a = sc.read_h5ad(f"{CKPT}/ckpt_14_clustered.h5ad")
print(f"loaded {a.shape}", flush=True)

SAMPLE_ORDER = ["N1-ABE","N2-ABE","N1-CBE","N2-CBE","N1-Control","N2-Control","ABE-Control","CBE-Control"]
# 语义独立 ⇒ 各用不同 family（同色=同义是读者的默认假设）
# pool = 二分类批次变量 ⇒ 用中性灰两档，把 8 个有彩色全部留给样本（pool 与 sample 同处 F3 一张图，不得共色）
C_POOL   = {"pool1_edited": PALETTE["grey_dark"], "pool2_control": PALETTE["grey_pale"]}
# phase：green_1(#DDF3DE) 白底对比度 1.17 不可见 ⇒ G2M 改用 neutral 浅灰
C_PHASE  = {"G1": PALETTE["green_2"], "S": PALETTE["green_3"], "G2M": PALETTE["neutral"]}
# 样本色：cool/warm 两 family 中白底可见的(去掉 #FFD700) + violet，恰 8 个
_cand = [c for c in (PALETTE_FAMILIES["cool"] + PALETTE_FAMILIES["warm"] + [PALETTE["violet"]])
         if contrast_vs_white(c) >= 1.5]
assert len(_cand) >= 8, f"可用样本色不足: {len(_cand)}"
C_SAMPLE = dict(zip(SAMPLE_ORDER, _cand[:8]))
assert len(set(C_SAMPLE.values())) == 8, "样本色有重复"
assert_no_palette_collision({"pool": C_POOL.values(), "phase": C_PHASE.values(), "sample": C_SAMPLE.values()})
_LEIDEN_BASE = [c for c in (PALETTE_FAMILIES["cool"] + PALETTE_FAMILIES["warm"] + [PALETTE["violet"]])
                if contrast_vs_white(c) >= 1.5] + [PALETTE["grey_soft"], PALETTE["grey_deep"], PALETTE["ink"]]
_LEIDEN_CATS = sorted(a.obs["leiden"].astype(str).unique(), key=lambda z: int(z) if z.isdigit() else z)
assert len(_LEIDEN_CATS) <= len(_LEIDEN_BASE), f"leiden {len(_LEIDEN_CATS)} 类 > 可用色 {len(_LEIDEN_BASE)}"
C_LEIDEN = dict(zip(_LEIDEN_CATS, _LEIDEN_BASE))
# 同图共存的组合：F3a = pool+sample；F3b/F4 = phase+leiden
assert_no_palette_collision({"phase": C_PHASE.values(), "leiden": C_LEIDEN.values()})
assert_visible_on_white(list(C_POOL.values()) + list(C_SAMPLE.values()) + list(C_PHASE.values())
                        + list(C_LEIDEN.values()), role="marker")
saved = []

# ---------- F1 QC 三联 ----------
fig, axes = create_subplots(1, 3, figsize=(13, 4))
for ax, (col, lab) in zip(axes, [("total_counts","UMI per cell"),
                                 ("n_genes_by_counts","Genes detected per cell"),
                                 ("pct_counts_mt","Mitochondrial % (11/13 MT genes only)")]):
    data = [a.obs.loc[a.obs.sample_id == s, col].values for s in SAMPLE_ORDER]
    bp = ax.boxplot(data, tick_labels=SAMPLE_ORDER, patch_artist=True, showfliers=False, widths=.6)
    for patch, s in zip(bp["boxes"], SAMPLE_ORDER):
        patch.set_facecolor(C_SAMPLE[s]); patch.set_edgecolor(PALETTE["ink"]); patch.set_linewidth(.8)
    for med in bp["medians"]: med.set_color(PALETTE["ink"]); med.set_linewidth(1.2)
    ax.set_ylabel(lab); ax.tick_params(axis="x", rotation=90)
    if col != "pct_counts_mt": ax.set_yscale("log")
axes[0].set_title("Per-sample QC after MAD filtering (293T, 10x Flex)", loc="left")
saved += finalize_figure(fig, f"{FIG}/F1_qc_per_sample.pdf", formats=["pdf"])

# ---------- F2 深度不均衡 ----------
col = pd.read_csv(f"{OUT}/16_pseudobulk_sample_coldata.csv", index_col=0).reindex(SAMPLE_ORDER)
fl  = pd.read_csv(f"{OUT}/12_per_sample_firstlook.csv").set_index("sample_id").reindex(SAMPLE_ORDER)
fig, axes = create_subplots(1, 2, figsize=(10, 4))
x = np.arange(len(SAMPLE_ORDER))
axes[0].bar(x, fl.reads_per_cell_vendor.values, color=[C_SAMPLE[s] for s in SAMPLE_ORDER],
            edgecolor=PALETTE["ink"], linewidth=.8)
axes[0].set_xticks(x); axes[0].set_xticklabels(SAMPLE_ORDER, rotation=90)
axes[0].set_ylabel("Vendor reads in cells / cell")
axes[0].set_title(f"Sequencing depth imbalance ({fl.reads_per_cell_vendor.max()/fl.reads_per_cell_vendor.min():.2f}x)",
                  loc="left")
axes[1].bar(x, col.median_umi.values, color=[C_SAMPLE[s] for s in SAMPLE_ORDER],
            edgecolor=PALETTE["ink"], linewidth=.8)
axes[1].set_xticks(x); axes[1].set_xticklabels(SAMPLE_ORDER, rotation=90)
axes[1].set_ylabel("Median UMI per cell (post-QC)")
axes[1].set_title("Not a biological difference", loc="left")
saved += finalize_figure(fig, f"{FIG}/F2_depth_imbalance.pdf", formats=["pdf"])

# ---------- F3 UMAP 四联（展示用；簇间距离不可定量解读）----------
um = a.obsm["X_umap"]
assert not np.isnan(um).any(), "UMAP 坐标含 NaN ⇒ 会静默丢点"
def _umap_panels(panels, out, suptitle):
    fig, axes = create_subplots(1, len(panels), figsize=(5.6 * len(panels), 5.2))
    for ax, (key, cmap, cats, title) in zip(axes, panels):
        drawn = 0
        for c in cats:
            m = (a.obs[key].astype(str) == c).values
            if m.sum() == 0: continue
            ax.scatter(um[m, 0], um[m, 1], s=.6, c=cmap[c], linewidths=0, rasterized=True, label=c)
            drawn += int(m.sum())
        assert drawn == a.n_obs, f"{key}: 画了 {drawn} 点，应为 {a.n_obs}（有类别被漏画）"
        ax.set_title(title, loc="left"); ax.set_xticks([]); ax.set_yticks([])
        ax.set_xlabel("UMAP1"); ax.set_ylabel("UMAP2")
        ax.legend(markerscale=12, frameon=False, fontsize=7, loc="best", ncol=1 if len(cats) <= 8 else 2)
    fig.suptitle(suptitle, fontsize=9, color=PALETTE["grey_dark"])
    return finalize_figure(fig, out, formats=["pdf"], pad=1.0)

_sup = "UMAP is for display only — inter-cluster distance carries no quantitative meaning"
saved += _umap_panels([("pool", C_POOL, list(C_POOL), "GEX pool (= batch, fully confounded with group)"),
                       ("sample_id", C_SAMPLE, SAMPLE_ORDER, "Sample")],
                      f"{FIG}/F3a_umap_pool_sample.pdf", _sup)
saved += _umap_panels([("phase", C_PHASE, list(C_PHASE), "Cell cycle phase"),
                       ("leiden", C_LEIDEN, _LEIDEN_CATS, "Leiden cluster")],
                      f"{FIG}/F3b_umap_phase_leiden.pdf", _sup)

# ---------- F4 组成：细胞周期 + 簇 ----------
ph = pd.read_csv(f"{OUT}/14_cellcycle_phase_pct.csv", index_col=0).reindex(SAMPLE_ORDER)
cs = pd.read_csv(f"{OUT}/14_cluster_by_sample_pct.csv", index_col=0).reindex(SAMPLE_ORDER)
fig, axes = create_subplots(1, 2, figsize=(12, 4.5))
bot = np.zeros(len(SAMPLE_ORDER))
for p in [c for c in ["G1","S","G2M"] if c in ph.columns]:
    axes[0].bar(x, ph[p].values, bottom=bot, color=C_PHASE[p], edgecolor=PALETTE["ink"],
                linewidth=.6, label=p)
    bot += ph[p].values
axes[0].set_xticks(x); axes[0].set_xticklabels(SAMPLE_ORDER, rotation=90)
axes[0].set_ylabel("% of cells"); axes[0].legend(frameon=False, fontsize=7)
axes[0].set_title("Cell cycle composition", loc="left")
bot = np.zeros(len(SAMPLE_ORDER))
for i, cl in enumerate(cs.columns):
    axes[1].bar(x, cs[cl].values, bottom=bot, color=C_LEIDEN[str(cl)],
                edgecolor=PALETTE["ink"], linewidth=.6, label=str(cl))
    bot += cs[cl].values
axes[1].set_xticks(x); axes[1].set_xticklabels(SAMPLE_ORDER, rotation=90)
axes[1].set_ylabel("% of cells")
axes[1].legend(frameon=False, fontsize=6, ncol=2, title="Leiden")
axes[1].set_title("Leiden cluster composition (not cell types — 293T is one cell line)", loc="left")
saved += finalize_figure(fig, f"{FIG}/F4_composition.pdf", formats=["pdf"])

# ---------- F5 pool 效应 vs pool 内效应（step 17 核心）----------
v = pd.read_csv(f"{OUT}/17_pool_effect_vs_within_pool_effect.csv")
bl = pd.read_csv(f"{OUT}/17_single_sample_pair_baseline.csv")
fig, (ax,) = create_subplots(1, 1, figsize=(8.5, 4.5))
lab = ["Confounded\ncross-pool\n(edited vs control)", "Identifiable\nwithin-pool\n(median)",
       "Single-sample pair\nbaseline (median)"]
cols = [PALETTE["red_strong"], PALETTE["blue_main"], PALETTE["grey_soft"]]
bars = ax.bar(np.arange(3), v.n_genes_abs_log2FC_gt_1.values, color=cols,
              edgecolor=PALETTE["ink"], linewidth=.9, width=.6)
for b, n in zip(bars, v.n_genes_abs_log2FC_gt_1.values):
    ax.text(b.get_x()+b.get_width()/2, b.get_height(), f"{int(n):,}", ha="center", va="bottom",
            fontsize=9, color=PALETTE["ink"])
ax.set_xticks(np.arange(3)); ax.set_xticklabels(lab, fontsize=8)
ax.set_ylabel("Genes with |log2FC| > 1")
_r = float(v.ratio_vs_within_pool_median.iloc[0])
ax.set_title(f"Confounded cross-pool contrast is {_r:.1f}× the within-pool median —\n"
             "group effect not separable from pool; no p-values reported", loc="left")
saved += finalize_figure(fig, f"{FIG}/F5_pool_vs_within_pool_effect.pdf", formats=["pdf"])

# ---------- F6 ambient rho ----------
rf = f"{OUT}/15_ambient_rho_percell.csv.gz"
if os.path.exists(rf):
    rho = pd.read_csv(rf)
    fig, (ax,) = create_subplots(1, 1, figsize=(8, 4))
    data = [rho.loc[rho.sample_id == s, "decontX_contamination"].values for s in SAMPLE_ORDER]
    bp = ax.boxplot(data, tick_labels=SAMPLE_ORDER, patch_artist=True, showfliers=False, widths=.6)
    for patch, s in zip(bp["boxes"], SAMPLE_ORDER):
        patch.set_facecolor(C_SAMPLE[s]); patch.set_edgecolor(PALETTE["ink"]); patch.set_linewidth(.8)
    for med in bp["medians"]: med.set_color(PALETTE["ink"]); med.set_linewidth(1.2)
    ax.set_ylabel("decontX contamination fraction")
    ax.tick_params(axis="x", rotation=90)
    ax.set_title("Ambient contamination (within-pool: includes other samples' transcripts)", loc="left")
    saved += finalize_figure(fig, f"{FIG}/F6_ambient_rho.pdf", formats=["pdf"])

pdfs = sorted(p for p in os.listdir(FIG) if p.endswith(".pdf"))
print(f"\n图件 {len(pdfs)} 个 PDF -> {FIG}")
for p in pdfs: print("   ", p)
assert not any(p.endswith(".png") for p in os.listdir(FIG)), "交付图件只要 PDF"
print("STEP19_PASS")
