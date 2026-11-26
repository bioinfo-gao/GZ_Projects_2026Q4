#!/usr/bin/env bash
# P47 round1 step 11: 0 算力静态预检。全部通过才允许进入 step 12。
# 判据来自 /scrna §5.1(静态扫描) 与 §2(已知静默坑)。
set -uo pipefail
PY=/mnt/ex_8T_SSD/conda_envs/sc_tut_main/bin/python
RS=/mnt/ex_8T_SSD/conda_envs/p47_r/bin/Rscript
FAIL=0
echo "##### 1. 解释器指纹（永远看 sys.executable，不看 CONDA_PREFIX）"
$PY -c "import sys;print('python :',sys.executable)" || FAIL=1
$RS -e 'cat("Rscript:", R.home(), "\n")' 2>/dev/null || FAIL=1

echo "##### 2. Python 侧会被调用的函数所在包（不止 library/import）"
$PY - <<'PY'
import sys, importlib
need = {"scanpy":["pp.scrublet","pp.highly_variable_genes","pp.neighbors","tl.leiden","tl.score_genes_cell_cycle","external"],
        "anndata":["read_h5ad","concat"], "h5py":[], "numpy":[], "scipy":["sparse"],
        "pandas":[], "skimage":[], "matplotlib":[], "harmonypy":[]}
bad=[]
for m in need:
    try:
        mod=importlib.import_module(m)
        print(f"  OK  {m:12s} {getattr(mod,'__version__','?')}")
    except Exception as e:
        print(f"  MISS {m:12s} {type(e).__name__}: {e}"); bad.append(m)
# scanpy 的具体函数存在性（版本间会搬家）
import scanpy as sc
for p in ["pp.scrublet","tl.score_genes_cell_cycle","pp.normalize_total","pp.log1p","tl.leiden","tl.umap","pp.pca"]:
    o=sc
    try:
        for part in p.split("."): o=getattr(o,part)
        print(f"  OK  sc.{p}")
    except AttributeError:
        print(f"  MISS sc.{p}"); bad.append("sc."+p)
# scrublet threshold=None 硬依赖 skimage
try:
    import skimage.filters; print("  OK  skimage.filters (scrublet threshold=None 依赖)")
except Exception as e:
    print("  MISS skimage.filters:",e); bad.append("skimage.filters")
sys.exit(1 if bad else 0)
PY
[ $? -ne 0 ] && FAIL=1

echo "##### 3. R 侧包 + scDblFinder x xgboost 直接调用断言（skill 指名判据）"
$RS -e 'for (p in c("scDblFinder","SoupX","celda","Seurat","SingleCellExperiment","Matrix","xgboost")) {
  ok <- requireNamespace(p, quietly=TRUE)
  cat(sprintf("  %s %-20s %s\n", ifelse(ok,"OK ","MISS"), p, ifelse(ok, as.character(packageVersion(p)), "")))
}' 2>&1 | grep -v "^$"
echo "  --- .xgbtrain 直接调用（tryCatch 会静默吞错，只能这样测）:"
$RS /home/gao/projects_2026Q4/47_Liu_Yuan_scRNA/docs/P47_scRNAseq_best_practices_evidence_1104/check7_scdblfinder_xgbtrain_direct.R 2>&1 | sed 's/^/  /'

echo "##### 4. 输入暂存完整性"
DST=/mnt/ex_8T_SSD/projects_2026Q4/47_Liu_Yuan_scRNA/staged_inputs_1125
N=$(find "$DST" -name "*.h5" 2>/dev/null | wc -l)
echo "  staged h5: $N (expect 16)"; [ "$N" -eq 16 ] || { echo "  FAIL 暂存未完成"; FAIL=1; }

echo; if [ $FAIL -eq 0 ]; then echo "PRECHECK_PASS"; else echo "PRECHECK_FAIL"; exit 1; fi
