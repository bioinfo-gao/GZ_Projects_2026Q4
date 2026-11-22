#!/usr/bin/env bash
# install-conda.sh —— 用 conda-forge + bioconda 二进制包安装本教程全部依赖
# 原机：miniforge base 环境 R 4.5.3 中验证。
# 本服务器（2026-11-21 适配）：装进独立前缀环境，不动 base 和其他项目的环境；
#   R 固定 4.5.3（与 _freeze/ 里的 sessionInfo 一致）；Quarto CLI 也装进该环境（服务器上没有系统 Quarto）。
# 用法：./install-conda.sh [环境前缀路径，默认 /mnt/ex_8T_SSD/conda_envs/p53_rbio]
# 之后：mamba activate /mnt/ex_8T_SSD/conda_envs/p53_rbio && quarto render
set -euo pipefail
PREFIX="${1:-/mnt/ex_8T_SSD/conda_envs/p53_rbio}"

PKGS=(
  r-base=4.5.3 quarto
  r-tidyverse r-lobstr r-bench r-r6 r-s7 r-rcpp r-renv r-targets r-tarchetypes r-qs2 r-here r-fs
  r-janitor r-glue r-conflicted r-withr r-matrixstats r-biocmanager r-remotes r-profvis r-sessioninfo
  r-data.table r-plyr r-dtplyr r-collapse r-arrow r-duckdb r-dbplyr
  r-ggrepel r-patchwork r-ggridges r-viridis r-scales r-pheatmap r-ggvenndiagram r-upsetr r-cowplot
  r-hexbin r-msigdbr r-ashr r-poiclaclu r-circlize r-ggbeeswarm
  r-seurat r-harmony r-r.utils r-hdf5r r-presto
  r-survival r-ggsurvfit r-survminer r-gtsummary r-broom.helpers
  bioconductor-summarizedexperiment bioconductor-singlecellexperiment bioconductor-genomicranges
  bioconductor-biostrings bioconductor-annotationdbi bioconductor-annotationhub bioconductor-experimenthub
  bioconductor-org.hs.eg.db bioconductor-org.mm.eg.db bioconductor-biocparallel bioconductor-hdf5array
  bioconductor-airway bioconductor-pasilla bioconductor-tximeta bioconductor-tximport bioconductor-deseq2
  bioconductor-apeglm bioconductor-edger bioconductor-limma bioconductor-vsn bioconductor-sva
  bioconductor-complexheatmap bioconductor-enhancedvolcano
  bioconductor-clusterprofiler bioconductor-enrichplot bioconductor-dose bioconductor-fgsea bioconductor-reactomepa
  bioconductor-scuttle bioconductor-scran bioconductor-scater bioconductor-bluster bioconductor-scdblfinder
  bioconductor-singler bioconductor-celldex bioconductor-glmgampoi bioconductor-muscat bioconductor-muscdata
  bioconductor-scrnaseq bioconductor-slingshot bioconductor-tscan bioconductor-zellkonverter
  bioconductor-dropletutils bioconductor-maftools bioconductor-variantannotation
  r-naniar r-visdat r-hgnchelper r-testthat r-clustree r-soupx r-mice r-checkmate r-timeroc r-cmprsk r-slider
  bioconductor-geoquery bioconductor-batchelor bioconductor-gsva bioconductor-rtracklayer
  bioconductor-ensdb.hsapiens.v86 bioconductor-txdb.hsapiens.ucsc.hg38.knowngene bioconductor-biomart
)

if [ -d "$PREFIX/conda-meta" ]; then
  mamba install -y -p "$PREFIX" -c conda-forge -c bioconda "${PKGS[@]}"
else
  mamba create  -y -p "$PREFIX" -c conda-forge -c bioconda "${PKGS[@]}"
fi

# 本服务器的环境变量（只在激活此环境时生效）：
# - BLAS/OpenMP 线程：conda R 的 OpenBLAS 默认吃满 64 核，超过服务器 56 线程上限；
#   必须在 R 启动前设置（.Renviron 读得太晚，OpenBLAS 已初始化）。
# - R_USER_CACHE_DIR：ExperimentHub / AnnotationHub / BiocFileCache / biomaRt / msigdbr 的缓存
#   默认在 ~/.cache/R，而 /home 所在盘已近满，改到 8T SSD。
# - XDG_CACHE_HOME：~/.cache 是指向 /Work_bio/gao/configs/.cache 的软链；Quarto 1.11.9 的 Deno KV
#   经软链路径打不开 sass.kv（"unable to open database file"，R 代码跑完后在 Sass 阶段失败）。
#   设成软链的真实路径 = 同一个目录，只是不经过软链。
CACHE_DIR=/mnt/ex_8T_SSD/R_user_cache
mkdir -p "$CACHE_DIR"
# mamba 2.x 不支持 `env config vars`，这里必须用 conda
conda env config vars set -p "$PREFIX" \
  OPENBLAS_NUM_THREADS=8 OMP_NUM_THREADS=8 R_USER_CACHE_DIR="$CACHE_DIR" \
  XDG_CACHE_HOME="$(readlink -f "$HOME/.cache")"

# 注意：不要用 conda 安装 r-quarto（R 包）—— 它会附带 conda 版 quarto 命令行，可能遮蔽系统 Quarto。
# 本服务器没有系统 Quarto，所以上面直接装了 conda-forge 的 quarto CLI 到本环境。
"$PREFIX/bin/Rscript" -e 'BiocManager::valid()' || true
