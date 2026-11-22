# install.R —— 一次性安装本教程用到的全部包
# 用法：在项目根目录执行  Rscript install.R
# 如果你的 R 来自 conda/miniforge，请改用 ./install-conda.sh（二进制包，免编译，见附录 A）
# 需要 R >= 4.4（推荐 4.5+）。Linux 用户请先装好系统依赖（见附录 A）。

options(repos = c(CRAN = "https://cloud.r-project.org"), timeout = 1800)

if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")

cran <- c(
  # 语言与工具
  "tidyverse", "rlang", "lobstr", "bench", "R6", "S7", "Rcpp", "renv", "targets",
  "tarchetypes", "qs2", "here", "fs", "janitor", "glue", "conflicted", "withr", "sessioninfo",
  # 数据处理生态（第 3 章）
  "data.table", "plyr", "dtplyr", "collapse", "arrow", "duckdb", "dbplyr",
  # 生存分析（第 14 章）
  "survival", "ggsurvfit", "survminer", "gtsummary", "broom.helpers",
  # 可视化
  "ggrepel", "patchwork", "ggridges", "viridis", "scales", "pheatmap",
  "ggVennDiagram", "UpSetR", "cowplot", "hexbin",
  # 单细胞（CRAN 侧）
  "Seurat", "harmony", "R.utils", "leidenbase", "PoiClaClu",
  # 详细版新增：缺失值、ID 清洗、聚类树、环境 RNA、插补、断言、生存
  "naniar", "visdat", "HGNChelper", "testthat", "clustree", "SoupX", "mice", "checkmate",
  "timeROC", "cmprsk", "slider", "httr2"
)

bioc <- c(
  # 核心基础设施
  "SummarizedExperiment", "SingleCellExperiment", "GenomicRanges", "IRanges",
  "Biostrings", "AnnotationDbi", "AnnotationHub", "ExperimentHub", "org.Hs.eg.db",
  "org.Mm.eg.db", "BiocParallel", "HDF5Array", "DelayedArray",
  # Bulk RNA-seq
  "airway", "pasilla", "tximeta", "tximport", "DESeq2", "apeglm", "ashr", "edgeR",
  "limma", "vsn", "sva", "ComplexHeatmap", "EnhancedVolcano",
  # 功能富集
  "clusterProfiler", "enrichplot", "DOSE", "fgsea", "msigdbr", "ReactomePA",
  # 单细胞（Bioconductor 侧）
  "scuttle", "scran", "scater", "bluster", "scDblFinder", "SingleR", "celldex",
  "glmGamPoi", "muscat", "muscData", "scRNAseq", "slingshot", "TSCAN",
  "zellkonverter", "DropletUtils",
  # 肿瘤基因组（第 14 章）
  "maftools", "VariantAnnotation",
  # 详细版新增：参考基因组与注释、GEO、GSVA、批次校正
  "GEOquery", "biomaRt", "rtracklayer", "GenomeInfoDb", "EnsDb.Hsapiens.v86",
  "TxDb.Hsapiens.UCSC.hg38.knownGene", "GSVA", "batchelor"
)

BiocManager::install(c(cran, bioc), update = FALSE, ask = FALSE, Ncpus = max(1, parallel::detectCores() - 1))

# GitHub 上的包（presto：快速 Wilcoxon；BPCells：可选，第 15 章大数据示例）
if (!requireNamespace("presto", quietly = TRUE)) remotes::install_github("immunogenomics/presto")
if (!requireNamespace("BPCells", quietly = TRUE)) try(remotes::install_github("bnprks/BPCells/r"))

missing <- setdiff(c(cran, bioc), rownames(installed.packages()))
if (length(missing)) {
  message("以下包未安装成功，请单独排查：\n  ", paste(missing, collapse = ", "))
} else {
  message("全部安装完成 ✔")
}
BiocManager::valid()
