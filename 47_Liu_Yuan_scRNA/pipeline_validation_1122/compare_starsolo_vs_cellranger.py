#!/usr/bin/env python3
"""STARsolo（本机）vs 10x 官方 Cell Ranger 3.0（pbmc_1k_v3）独立互证（2026-11-22）
注释不同（CR 3.0 = GRCh38-3.0.0/Ensembl 93；本机 = GENCODE v45 + 10x-2024A 规则过滤）⇒ 预期高度相关、不完全相同。
CR 3.0 只计外显子 ⇒ 对比 STARsolo 的 Gene（不是 GeneFull_Ex50pAS）。
"""
import sys, numpy as np, scanpy as sc, scipy.stats as ss
solo, cr = sys.argv[1], sys.argv[2]   # Solo.out/Gene/filtered  ,  filtered_feature_bc_matrix
a = sc.read_10x_mtx(solo, var_names="gene_ids"); b = sc.read_10x_mtx(cr, var_names="gene_ids")
b.obs_names = [x.split("-")[0] for x in b.obs_names]; a.obs_names = [x.split("-")[0] for x in a.obs_names]
strip = lambda v: [x.split(".")[0] for x in v]
a.var_names = strip(a.var_names); b.var_names = strip(b.var_names)
a.var_names_make_unique(); b.var_names_make_unique()
ca, cb = set(a.obs_names), set(b.obs_names); both = sorted(ca & cb)
print(f"cells: STARsolo={len(ca)}  CellRanger={len(cb)}  shared={len(both)}  "
      f"shared/CR={len(both)/len(cb):.1%}  shared/STARsolo={len(both)/len(ca):.1%}")
ua = np.asarray(a[both].X.sum(1)).ravel(); ub = np.asarray(b[both].X.sum(1)).ravel()
print(f"per-cell UMI (shared cells): median STARsolo={np.median(ua):.0f}  CR={np.median(ub):.0f}  "
      f"ratio={np.median(ua/ub):.3f}  spearman={ss.spearmanr(ua, ub)[0]:.4f}")
ga = np.asarray((a[both].X > 0).sum(1)).ravel(); gb = np.asarray((b[both].X > 0).sum(1)).ravel()
print(f"per-cell genes: median STARsolo={np.median(ga):.0f}  CR={np.median(gb):.0f}")
g = sorted(set(a.var_names) & set(b.var_names))
ta = np.asarray(a[both, g].X.sum(0)).ravel(); tb = np.asarray(b[both, g].X.sum(0)).ravel()
m = (ta + tb) > 0
print(f"genes: shared IDs={len(g)}  pearson(log1p total UMI)={np.corrcoef(np.log1p(ta[m]), np.log1p(tb[m]))[0,1]:.4f}")
# 每细胞基因表达谱一致性（共有细胞 × 共有基因，log1p CPM 的 Pearson 中位数）
def cpm(x): x = np.asarray(x.todense()) if hasattr(x, "todense") else x; return np.log1p(x / x.sum(1, keepdims=True) * 1e4)
xa, xb = cpm(a[both, g].X), cpm(b[both, g].X)
r = [np.corrcoef(xa[i], xb[i])[0, 1] for i in range(len(both))]
print(f"per-cell profile pearson: median={np.median(r):.4f}  5th pct={np.percentile(r, 5):.4f}")
