#!/usr/bin/env python
"""P47 round1 step 12b：对【已落盘的 ckpt_12】重跑体检断言（修正期望值后）。
为什么不重跑 step 12：step 12 的 checkpoint 在断言失败之前就已写完且内容正确，
失败原因是我把"矩阵 MT 基因数"的期望值写成了面板的 12（实际 11，MT-CO1 included=FALSE）。
验证落盘产物本身，比验证代码路径更强；也避免重写 17.5 GB。
"""
import os, sys, json
for v in ("OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS","NUMEXPR_NUM_THREADS"):
    os.environ[v]="16"
import numpy as np, pandas as pd, scanpy as sc
CKPT="/mnt/ex_8T_SSD/projects_2026Q4/47_Liu_Yuan_scRNA/ckpt_1125"
STAGE="/mnt/ex_8T_SSD/projects_2026Q4/47_Liu_Yuan_scRNA/staged_inputs_1125"
OUT="/home/gao/projects_2026Q4/47_Liu_Yuan_scRNA/output_results/round1_descriptive_1125"
print("python :", sys.executable, flush=True)
a = sc.read_h5ad(f"{CKPT}/ckpt_12_raw_concat.h5ad")
print("loaded", a.shape, flush=True)
EXPECT={"N1-ABE":15245,"N2-ABE":19027,"N1-CBE":20706,"N2-CBE":15726,
        "N1-Control":20780,"N2-Control":20053,"CBE-Control":12480,"ABE-Control":12027}
ps=pd.read_csv(f"{STAGE}/N1-ABE-N2-ABE-N1-CBE-N2-CBE/probe_set.csv",comment="#")
ps["g"]=ps.probe_id.str.split("|").str[1]
inc=ps[ps.included.astype(str).str.upper()=="TRUE"]
n_panel, n_inc = ps.gene_id.nunique(), inc.gene_id.nunique()
mt_mat=sorted(g for g in a.var_names if g.startswith("MT-"))
fails=[]
def ck(name, ok, got, want):
    print(f"  {'OK ' if ok else 'FAIL'} {name}: got={got} want={want}")
    if not ok: fails.append(name)
ck("total cells", a.n_obs==sum(EXPECT.values()), a.n_obs, sum(EXPECT.values()))
ck("features == genes with >=1 included probe", a.n_vars==n_inc, a.n_vars, n_inc)
ck("992 genes dropped (all probes FALSE)", n_panel-n_inc==992, n_panel-n_inc, 992)
ck("MT genes in matrix", len(mt_mat)==11, len(mt_mat), "11 (= panel 12 - MT-CO1)")
ck("MT-CO1 absent from matrix", "MT-CO1" not in a.var_names, "MT-CO1" in a.var_names, "absent")
ck("MT-ATP8 absent from panel", "MT-ATP8" not in set(ps.g), "MT-ATP8" in set(ps.g), "absent")
ck("obs_names unique", a.obs_names.is_unique, a.obs_names.is_unique, True)
ck("per-sample cells match vendor", dict(a.obs.sample_id.value_counts())==EXPECT, "see table", "match")
ck("barcode core length 26", sorted({len(b.split('_')[-1].split('-')[0]) for b in a.obs_names[:5000]})==[26],
   sorted({len(b.split('_')[-1].split('-')[0]) for b in a.obs_names[:5000]}), [26])
ck("layers['counts'] present & integral", "counts" in a.layers and
   bool(np.all(a.layers["counts"].data==np.round(a.layers["counts"].data))), "integral", "integral")
ck("pct_counts_mt > 0 for all samples",
   bool(a.obs.groupby("sample_id",observed=True).pct_counts_mt.mean().gt(0).all()), "all>0", "all>0")
nnz=int(a.layers["counts"].nnz)
print(f"\n### 稠密度：nnz={nnz:,}  = {100*nnz/(a.n_obs*a.n_vars):.1f}% of matrix"
      f"  ⇒ 平均每细胞检出 {nnz/a.n_obs:,.0f} 基因 / {a.n_vars}")
mtp = a.obs.groupby("sample_id",observed=True).pct_counts_mt.describe()[["mean","50%","max"]].round(3)
print("\n### pct_counts_mt（只基于 11/13 个 MT 基因，不可与常规 scRNA 经验阈值相比）:")
print(mtp.to_string())
json.dump({"n_panel_genes":int(n_panel),"n_features_in_matrix":int(n_inc),
           "n_genes_all_probes_excluded":int(n_panel-n_inc),
           "mt_genes_in_matrix":mt_mat,"mt_genes_quantified_of_human13":len(mt_mat),
           "mt_absent_reason":{"MT-CO1":"probe included=FALSE","MT-ATP8":"not in panel"},
           "nnz":nnz,"density_pct":round(100*nnz/(a.n_obs*a.n_vars),3),
           "mean_genes_detected_per_cell":round(nnz/a.n_obs,1),
           "pct_mt_not_comparable_to_standard_scRNA":True},
          open(f"{OUT}/12b_matrix_structure_facts.json","w"), indent=2)
mtp.to_csv(f"{OUT}/12b_pct_mt_per_sample.csv")
if fails: print("\n#### 失败 ####\n"+"\n".join(fails)); sys.exit(1)
print("\nSTEP12B_PASS")
