#!/usr/bin/env python
"""P47 round1 step 12：读入 8 个样本并做第一次体检（Flex 特异）。

上游：staged_inputs_1125/（step 10）  产出：ckpt_12_raw_concat.h5ad + 体检审计表
判据来源：/scrna §3.1，但【MT 基因数那条必须改】——人类常规期望 13，
本批 Flex 面板实测只有 12（缺 MT-ATP8），照搬 13 会误判。
"""
import os, sys, json
for v in ("OMP_NUM_THREADS","OPENBLAS_NUM_THREADS","MKL_NUM_THREADS","NUMEXPR_NUM_THREADS"):
    os.environ[v] = "16"          # P46 DELLY 在跑，全机 28 核硬顶需共享
import numpy as np, pandas as pd, scanpy as sc, anndata as ad

STAGE = "/mnt/ex_8T_SSD/projects_2026Q4/47_Liu_Yuan_scRNA/staged_inputs_1125"
OUT   = "/home/gao/projects_2026Q4/47_Liu_Yuan_scRNA/output_results/round1_descriptive_1125"
CKPT  = "/mnt/ex_8T_SSD/projects_2026Q4/47_Liu_Yuan_scRNA/ckpt_1125"
for d in (OUT, CKPT): os.makedirs(d, exist_ok=True)
print("python :", sys.executable)

POOLS = {
    "pool1_edited":  "N1-ABE-N2-ABE-N1-CBE-N2-CBE",
    "pool2_control": "N1-Control-N2-Control-CBE-Control-ABE-Control",
}
# 公司 QC 的 cell 数，用作读入正确性的独立互证（不是我方重算值）
EXPECT_CELLS = {"N1-ABE":15245,"N2-ABE":19027,"N1-CBE":20706,"N2-CBE":15726,
                "N1-Control":20780,"N2-Control":20053,"CBE-Control":12480,"ABE-Control":12027}
GROUP = {s:("edited" if s in ("N1-ABE","N2-ABE","N1-CBE","N2-CBE") else "control") for s in EXPECT_CELLS}
EDITOR = {"N1-ABE":"ABE","N2-ABE":"ABE","N1-CBE":"CBE","N2-CBE":"CBE",
          "N1-Control":"none_gRNA_ctrl","N2-Control":"none_gRNA_ctrl",
          "CBE-Control":"CBE_editor_only","ABE-Control":"ABE_editor_only"}
GUIDE = {"N1-ABE":"N1","N2-ABE":"N2","N1-CBE":"N1","N2-CBE":"N2",
         "N1-Control":"N1","N2-Control":"N2","CBE-Control":"none","ABE-Control":"none"}

fails, rows = [], []
def check(name, ok, got, want):
    (print if ok else print)(f"  {'OK ' if ok else 'FAIL'} {name}: got={got} want={want}")
    if not ok: fails.append(f"{name}: got={got} want={want}")

# ---------- 面板事实：MT 基因到底有几个 ----------
ps = pd.read_csv(f"{STAGE}/{POOLS['pool1_edited']}/probe_set.csv", comment="#")
ps["gene_name"] = ps["probe_id"].str.split("|").str[1]
mt_panel = sorted(ps.loc[ps.gene_name.str.startswith("MT-"), "gene_name"].unique())
print(f"\n### 面板 MT 基因：{len(mt_panel)} 个 -> {mt_panel}")
check("probe panel MT gene count", len(mt_panel) == 12, len(mt_panel), "12 (Flex 面板事实，非人类常规 13)")
HUMAN13 = {"MT-ND1","MT-ND2","MT-CO1","MT-CO2","MT-ATP8","MT-ATP6","MT-CO3",
           "MT-ND3","MT-ND4L","MT-ND4","MT-ND5","MT-ND6","MT-CYB"}
print(f"### 相对人类常规 13 个缺少：{sorted(HUMAN13 - set(mt_panel))}")

adatas = {}
for pool, lib in POOLS.items():
    base = f"{STAGE}/{lib}/per_sample_outs"
    for smp in sorted(os.listdir(base)):
        h5 = f"{base}/{smp}/sample_filtered_feature_bc_matrix.h5"
        a = sc.read_10x_h5(h5)
        n_dup_before = int(a.var_names.duplicated().sum())
        a.var_names_make_unique()                      # 原地修改，返回 None
        bclen = sorted({len(b.split("-")[0]) for b in a.obs_names})
        ftypes = sorted(set(a.var.get("feature_types", pd.Series(["<absent>"])).astype(str)))
        mt_in_mat = sorted([g for g in a.var_names if g.startswith("MT-")])
        print(f"\n### {smp}  ({pool})  shape={a.shape}")
        check(f"{smp} n_cells vs vendor QC", a.n_obs == EXPECT_CELLS[smp], a.n_obs, EXPECT_CELLS[smp])
        check(f"{smp} barcode core length", bclen == [26], bclen, "[26] (16 GEM + 10 probe BC)")
        check(f"{smp} feature_types unique", ftypes == ["Gene Expression"], ftypes, "['Gene Expression']")
        # ⚠ 矩阵里是 11 不是 12：MT-CO1 的唯一探针 included=FALSE，被排除出 filtered probe set。
        # 面板 12 个 MT 基因 - MT-CO1 = 11；再算上面板本就没有 MT-ATP8，
        # ⇒ 人类 13 个线粒体蛋白编码基因里只有 11 个被定量。
        # 后果：本数据的 pct_counts_mt 不可与常规 scRNA 的 mito% 经验阈值相比。
        check(f"{smp} MT genes in matrix", len(mt_in_mat) == 11, len(mt_in_mat),
              "11 (= 面板 12 - MT-CO1(included=FALSE))")
        a.var["mt"] = a.var_names.str.startswith("MT-")
        sc.pp.calculate_qc_metrics(a, qc_vars=["mt"], inplace=True, percent_top=None, log1p=False)
        check(f"{smp} pct_counts_mt not all zero", float(a.obs.pct_counts_mt.sum()) > 0,
              round(float(a.obs.pct_counts_mt.mean()), 3), "> 0 (全 0 = 前缀写错的典型静默坑)")
        a.obs["sample_id"], a.obs["pool"] = smp, pool
        a.obs["group"], a.obs["editor"], a.obs["guide"] = GROUP[smp], EDITOR[smp], GUIDE[smp]
        a.obs_names = [f"{smp}_{b}" for b in a.obs_names]      # 显式前缀，不依赖 concat 后缀
        rows.append(dict(sample_id=smp, pool=pool, group=GROUP[smp], editor=EDITOR[smp], guide=GUIDE[smp],
                         n_cells=a.n_obs, n_genes_detected=int((a.var.n_cells_by_counts > 0).sum()),
                         n_features_in_h5=a.n_vars, dup_var_names_before_unique=n_dup_before,
                         median_umi=float(np.median(a.obs.total_counts)),
                         median_genes=float(np.median(a.obs.n_genes_by_counts)),
                         total_umi=float(a.obs.total_counts.sum()),
                         mean_pct_mt=float(a.obs.pct_counts_mt.mean()),
                         max_pct_mt=float(a.obs.pct_counts_mt.max())))
        adatas[smp] = a

# ---------- 合并：按 gene_ids 不按 symbol ----------
gid_sets = {s: set(a.var["gene_ids"]) for s, a in adatas.items()}
common = set.intersection(*gid_sets.values())
check("all samples share identical gene_ids set", all(len(v) == len(common) for v in gid_sets.values()),
      sorted({len(v) for v in gid_sets.values()}), f"all == {len(common)}")
adata = ad.concat(adatas, label="sample_key", index_unique=None, join="outer")
adata.var = adatas["N1-ABE"].var[["gene_ids","feature_types","mt"]].copy()
check("concat total cells", adata.n_obs == sum(EXPECT_CELLS.values()), adata.n_obs, sum(EXPECT_CELLS.values()))
check("concat obs_names unique", adata.obs_names.is_unique, adata.obs_names.is_unique, True)
adata.layers["counts"] = adata.X.copy()            # 归一化之前显式备份
check("layers['counts'] is integer", bool(np.all(adata.layers["counts"].data == np.round(adata.layers["counts"].data))),
      "integral", "integral")

df = pd.DataFrame(rows).sort_values(["pool","sample_id"])
df["reads_in_cells_vendor"] = df.sample_id.map({"N1-ABE":1358931130,"N2-ABE":1859024169,"N1-CBE":1941921733,
    "N2-CBE":1693758770,"N1-Control":2222181629,"N2-Control":2192306043,"CBE-Control":746751973,"ABE-Control":949000266})
df["reads_per_cell_vendor"] = (df.reads_in_cells_vendor / df.n_cells).round(0)
df.to_csv(f"{OUT}/12_per_sample_firstlook.csv", index=False)

depth_ratio = df.reads_in_cells_vendor.max() / df.reads_in_cells_vendor.min()
umi_ratio   = df.median_umi.max() / df.median_umi.min()
print(f"\n### 深度不均衡：vendor reads-in-cells 极差 {depth_ratio:.2f}× ；median UMI/cell 极差 {umi_ratio:.2f}×")
pd.DataFrame([dict(metric="vendor_reads_in_cells_max_over_min", value=round(float(depth_ratio),3)),
              dict(metric="median_umi_max_over_min",            value=round(float(umi_ratio),3)),
              dict(metric="panel_mt_gene_count",                value=len(mt_panel)),
              dict(metric="total_cells",                        value=int(adata.n_obs)),
              dict(metric="n_features",                         value=int(adata.n_vars))
              ]).to_csv(f"{OUT}/12_imbalance_and_panel_facts.csv", index=False)
json.dump({"panel_mt_genes": mt_panel, "missing_vs_human13": sorted(HUMAN13 - set(mt_panel)),
           "depth_ratio_vendor_reads": float(depth_ratio), "umi_ratio_median": float(umi_ratio)},
          open(f"{OUT}/12_panel_and_depth_facts.json","w"), indent=2)

adata.write(f"{CKPT}/ckpt_12_raw_concat.h5ad")
print(f"\nwrote {CKPT}/ckpt_12_raw_concat.h5ad  shape={adata.shape}")
print(df[["sample_id","pool","group","n_cells","median_umi","median_genes","mean_pct_mt","reads_per_cell_vendor"]].to_string(index=False))
if fails:
    print("\n".join(["", "#### 断言失败 ####"] + fails)); sys.exit(1)
print("\nSTEP12_PASS")
