"""pbmc_1k_v3（10x 3' v3，~54k reads/cell）chrM 每 100 bp 网格的分子层面覆盖（2026-11-22）"""
import pandas as pd, numpy as np
s = pd.read_csv("pbmc_chrM_grid100.site_audit.tsv", sep="\t")
m = pd.read_csv("pbmc_chrM_grid100.molecules.tsv", sep="\t")
ncell = sum(1 for _ in open("cells.tsv"))
s["sec_frac"] = s.rec_secondary / s.rec_all.replace(0, np.nan)
s["cell_frac"] = s.cells_with_molecule / ncell
pc = m.groupby(["label", "CB"]).size().rename("mol").reset_index()
med_pc = pc.groupby("label").mol.median()
s["mol_per_covered_cell_median"] = s.label.map(med_pc).fillna(0)
q = lambda x: f"min={x.min():.3g} p25={x.quantile(.25):.3g} median={x.median():.3g} p75={x.quantile(.75):.3g} max={x.max():.3g}"
print(f"QC cells={ncell}  grid sites={len(s)}")
print("molecules per site      :", q(s.molecules))
print("fraction of cells with >=1 molecule per site:", q(s.cell_frac))
print("molecules per covered cell (per-site median):", q(s.mol_per_covered_cell_median))
print("secondary record fraction per site:", q(s.sec_frac.fillna(0)))
for t in (0.1, 0.5, 0.9): print(f"sites where >= {t:.0%} of cells have >=1 molecule: {(s.cell_frac >= t).sum()}/{len(s)}")
# 每个细胞：166 个网格位点里有几个位点拿到 >=1 / >=5 个分子
cs = pc.groupby("CB").agg(sites1=("mol", lambda x: (x >= 1).sum()), sites5=("mol", lambda x: (x >= 5).sum()), mol_mean=("mol", "sum"))
cs["mol_mean"] = cs.mol_mean / len(s)
cs = cs.reindex(pd.read_csv("cells.tsv", header=None)[0]).fillna(0)
print("per cell: grid sites with >=1 molecule:", q(cs.sites1), f" (of {len(s)})")
print("per cell: grid sites with >=5 molecules:", q(cs.sites5))
print("per cell: mean molecules per grid site (= molecule-level depth):", q(cs.mol_mean))
s.sort_values("sec_frac", ascending=False).head(8)[["label", "rec_all", "rec_secondary", "sec_frac", "excl_mapq<20"]].to_csv("top_secondary_sites.tsv", sep="\t", index=False)
print("top secondary-fraction sites ->", open("top_secondary_sites.tsv").read())
