"""maegatk（bwa 重比对到内置 rCRS、UMI 共识）vs 5_target_site_audit.py（STAR 比对、primary MAPQ>=20、BQ>=20、(CB,UB) 分子）
两种实现、两种过滤口径 ⇒ 只期望高度相关，不期望数值相同。这只验证算法层面的一致性，不是独立生物学验证。"""
import glob, pandas as pd, numpy as np, scipy.stats as ss
cov_f = glob.glob("maegatk_out/final/*coverage.txt.gz")
assert len(cov_f) == 1, cov_f
cov = pd.read_csv(cov_f[0], header=None, names=["pos", "CB", "depth"])
s = pd.read_csv("pbmc_chrM_grid100.site_audit.tsv", sep="\t")
m = pd.read_csv("pbmc_chrM_grid100.molecules.tsv", sep="\t")
cells = set(pd.read_csv("cells.tsv", header=None)[0])
print(f"maegatk: positions={cov.pos.nunique()} cells={cov.CB.nunique()} (QC cells={len(cells)}; in QC set={len(set(cov.CB) & cells)})")
g = cov[cov.pos.isin(s.pos)]
site_mae = g.groupby("pos").depth.sum()
site_aud = s.set_index("pos").molecules
j = pd.concat([site_mae.rename("maegatk"), site_aud.rename("audit")], axis=1).fillna(0)
print(f"per-site total (166 grid sites): spearman={ss.spearmanr(j.maegatk, j.audit)[0]:.4f}  "
      f"median ratio maegatk/audit={np.median(j.maegatk[j.audit > 0] / j.audit[j.audit > 0]):.3f}")
pc_aud = m.groupby(["pos", "CB"]).size().rename("audit")
pc_mae = g.set_index(["pos", "CB"]).depth.rename("maegatk")
k = pd.concat([pc_mae, pc_aud], axis=1).fillna(0)
print(f"per site x cell ({len(k):,} pairs with coverage in either): spearman={ss.spearmanr(k.maegatk, k.audit)[0]:.4f}  "
      f"both>0: {((k.maegatk > 0) & (k.audit > 0)).mean():.1%}  only maegatk: {((k.maegatk > 0) & (k.audit == 0)).mean():.1%}  only audit: {((k.maegatk == 0) & (k.audit > 0)).mean():.1%}")
# 分歧最大的位点（看是否集中在 NUMT 高 secondary 区）
s["sec_frac"] = s.rec_secondary / s.rec_all.replace(0, np.nan)
j["log2_ratio"] = np.log2((j.maegatk + 1) / (j.audit + 1))
j = j.join(s.set_index("pos")[["sec_frac", "excl_mapq<20"]])
print("sites with largest maegatk/audit disagreement:")
print(j.reindex(j.log2_ratio.abs().sort_values(ascending=False).index).head(8).to_string())
print(f"spearman(|log2 ratio|, secondary fraction) = {ss.spearmanr(j.log2_ratio.abs(), j.sec_frac.fillna(0))[0]:.3f}")
