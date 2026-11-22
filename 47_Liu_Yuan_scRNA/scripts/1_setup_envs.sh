#!/usr/bin/env bash
# 1_setup_envs.sh —— P47 专用的两个 conda 环境（2026-11-21 建）
# 不动 sc_tut_r / sc_tut_main（已验证的教程环境）。表达主线仍用 sc_tut_main（scanpy 1.11.5）。
#   p47_r      : ambient RNA（SoupX / DecontX）、doublet（scDblFinder）、pseudobulk DE、读 h5ad
#   p47_allele : 靶位点按 barcode/UMI 的 allele 计数（cellsnp-lite / maegatk）
# 版本取舍：scDblFinder 1.26.9（Bioc 3.23）与 celda（R 4.5 只有 Bioc 3.22 build）装不进同一环境；
#   上游 1.27.9 修复的是 xgboost >=3.4 的静默回退问题 ⇒ 这里保留 1.24.10 并钉 r-xgboost<3.4（解出 3.3.1）。
#   正式运行仍要跑 scDblFinder 阳性对照、确认没有回退（见 docs/P47_scRNAseq_best_practices_github_survey_1104.md §5.3）。
# 用法：bash scripts/1_setup_envs.sh            # 新建（已存在则跳过）
#       bash scripts/1_setup_envs.sh --verify-only
set -euo pipefail
ROOT=/mnt/ex_8T_SSD/conda_envs
R_ENV=$ROOT/p47_r
A_ENV=$ROOT/p47_allele
CH=(-c conda-forge -c bioconda)

R_PKGS=(r-base=4.5.3 'r-xgboost<3.4' bioconductor-scdblfinder r-soupx bioconductor-celda
        bioconductor-singlecellexperiment bioconductor-scuttle bioconductor-scater
        bioconductor-deseq2 bioconductor-edger bioconductor-zellkonverter bioconductor-dropletutils
        r-hdf5r r-tidyverse)
A_PKGS=(cellsnp-lite=1.2.3 maegatk=0.2.0 samtools pysam)

if [ "${1:-}" != "--verify-only" ]; then
  [ -d "$R_ENV/conda-meta" ] || mamba create -y -p "$R_ENV" "${CH[@]}" "${R_PKGS[@]}"
  [ -d "$A_ENV/conda-meta" ] || mamba create -y -p "$A_ENV" "${CH[@]}" "${A_PKGS[@]}"
  # 本服务器：BLAS/OMP 线程上限 + R 缓存移出近满的 /home（mamba 2.x 无 env config vars，用 conda）
  conda env config vars set -p "$R_ENV" OPENBLAS_NUM_THREADS=8 OMP_NUM_THREADS=8 \
    R_USER_CACHE_DIR=/mnt/ex_8T_SSD/R_user_cache

  # 2026-11-22 增补：亚群校准（ClusterDE / count splitting），见 docs/P47_readiness_gaps_1122.md §3
  # conda 上没有 ClusterDE/scDesign3/countsplit ⇒ 依赖与编译器走 conda，这 5 个包源码装；GitHub 钉 SHA。
  # 编译器必须是 *_linux-64（R Makeconf 用 x86_64-conda-linux-gnu-cc）。
  mamba install -y -p "$R_ENV" "${CH[@]}" gcc_linux-64 gxx_linux-64 gfortran_linux-64 make \
    r-remotes r-rcpp r-rcpparmadillo r-mgcv r-gamlss r-gamlss.dist r-mclust r-mvtnorm r-pbmcapply \
    r-umap r-irlba r-igraph r-rvinecopulib r-mvnfast r-coop r-sparsemvn r-fitdistrplus r-paireddata \
    r-zigg r-gson bioconductor-cellmixs r-checkmate r-lsa r-ade4   # 后 3 个：bettermc / FRmatch 的依赖（首轮漏，1122 补）
  (eval "$(mamba shell hook --shell bash)"; mamba activate "$R_ENV"
   MAKEFLAGS=-j8 Rscript -e '
    options(repos = c(CRAN = "https://cloud.r-project.org"))
    if (!requireNamespace("countsplit", quietly = TRUE)) install.packages("countsplit")
    gh <- c("gfkse/bettermc@f11b9bc654d5", "NLM-DIR/FRmatch@edcb0253fbf8",
            "SONGDONGYUAN1994/scDesign3@dec8acb2f54c", "SONGDONGYUAN1994/ClusterDE@93d8c22599c9")
    for (r in gh) if (!requireNamespace(sub("@.*", "", sub(".*/", "", r)), quietly = TRUE))
      remotes::install_github(r, dependencies = FALSE, upgrade = "never")')
fi

echo "== verify $R_ENV"
"$R_ENV/bin/Rscript" -e '
pk <- c("scDblFinder","xgboost","SoupX","celda","SingleCellExperiment","scuttle","scater",
        "DESeq2","edgeR","zellkonverter","DropletUtils","hdf5r","tidyverse",
        "countsplit","scDesign3","ClusterDE")
v <- vapply(pk, function(p) if (requireNamespace(p, quietly = TRUE)) as.character(packageVersion(p)) else NA_character_, "")
print(v); stopifnot(!anyNA(v), package_version(v[["xgboost"]]) < "3.4")'
echo "== verify $A_ENV"
"$A_ENV/bin/cellsnp-lite" --version
"$A_ENV/bin/maegatk" --version
"$A_ENV/bin/samtools" --version | head -1
"$A_ENV/bin/python" -c "import pysam; print('pysam', pysam.__version__)"
echo VERIFY_OK
