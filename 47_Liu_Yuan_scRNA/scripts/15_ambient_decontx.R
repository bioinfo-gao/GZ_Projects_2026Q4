#!/usr/bin/env Rscript
# P47 round1 step 15：ambient RNA 污染率估计（decontX）。
#
# 为什么放在聚类之后：/scrna §3.2 —— 估 rho 要先定义"不该表达该基因的细胞"，
#   这个判据必须独立于被测基因本身。拿基因自己的分位数圈背景 = 循环论证，
#   结果按构造恒为 0（proj22 实测 6 个示踪基因 5 个返回 0.000）。decontX 用簇标签满足这一点。
# 为什么用 decontX 而不是 SoupX：SoupX 的 autoEstCont 估 soup 谱只用总 UMI 落在
#   soupRange=c(0,100) 的液滴；本批 Flex 的 raw 矩阵液滴 UMI 分布未必覆盖该窗口，
#   没有这样的液滴时 alloc() 直接报 "missing value where TRUE/FALSE needed"。
#   decontX 可直接吃 background=raw 矩阵，更稳。两者本机均已功能冒烟通过(2026-11-21)。
# ⚠ Flex 特异：同 pool 内 4 个样本共用 GEM，GEM 重叠实测 11.4–17.4%
#   ⇒ 这里的"ambient"同时包含【同 pool 其他样本的转录本】，不是单样本 soup。
suppressPackageStartupMessages({library(celda); library(SingleCellExperiment); library(Seurat); library(Matrix)})
set.seed(20261125)

STAGE <- "/mnt/ex_8T_SSD/projects_2026Q4/47_Liu_Yuan_scRNA/staged_inputs_1125"
OUT   <- "/home/gao/projects_2026Q4/47_Liu_Yuan_scRNA/output_results/round1_descriptive_1125"
cat("celda", as.character(packageVersion("celda")), "\n")

meta <- read.csv(gzfile(file.path(OUT, "14_percell_metadata.csv.gz")), row.names = 1)
stopifnot(nrow(meta) > 0, "leiden" %in% colnames(meta), "sample_id" %in% colnames(meta))
cat("per-cell metadata rows:", nrow(meta), " clusters:", length(unique(meta$leiden)), "\n")

libs <- c("N1-ABE-N2-ABE-N1-CBE-N2-CBE", "N1-Control-N2-Control-CBE-Control-ABE-Control")
rows <- list(); percell <- list()
for (lib in libs) {
  base <- file.path(STAGE, lib, "per_sample_outs")
  for (smp in sort(list.dirs(base, recursive = FALSE, full.names = FALSE))) {
    mf <- Read10X_h5(file.path(base, smp, "sample_filtered_feature_bc_matrix.h5"))
    mr <- Read10X_h5(file.path(base, smp, "sample_raw_feature_bc_matrix.h5"))
    # 只保留通过 QC 的细胞，并按名字对齐簇标签（绝不按位置赋值）
    sm <- meta[meta$sample_id == smp, , drop = FALSE]
    bc <- sub(paste0("^", smp, "_"), "", rownames(sm))
    keep <- bc %in% colnames(mf)
    stopifnot(all(keep))
    mf <- mf[, bc, drop = FALSE]
    z <- as.integer(factor(sm$leiden))
    names(z) <- bc
    stopifnot(identical(colnames(mf), names(z)), !any(is.na(z)))
    # 分母非零断言：背景液滴必须真的存在，否则 rho 会退化成一个长得像结论的数
    bg_only <- setdiff(colnames(mr), colnames(mf))
    bg_counts <- sum(mr[, head(bg_only, 20000), drop = FALSE])
    stopifnot(length(bg_only) > 1000, bg_counts > 0)

    sce <- SingleCellExperiment(list(counts = mf))
    t0 <- Sys.time()
    if (length(unique(z)) >= 2) {
      sce <- decontX(sce, z = z, background = SingleCellExperiment(list(counts = mr[, bg_only])), verbose = FALSE)
      zmode <- sprintf("leiden(%d levels)", length(unique(z)))
    } else {
      sce <- decontX(sce, background = SingleCellExperiment(list(counts = mr[, bg_only])), verbose = FALSE)
      zmode <- "decontX-internal(sample had <2 leiden levels)"
    }
    el <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1)
    rho <- sce$decontX_contamination
    stopifnot(length(unique(round(rho, 6))) > 100)   # 常数 rho = 估计退化
    cat(sprintf("  %-12s cells=%6d  rho median=%.4f mean=%.4f p90=%.4f  bg droplets=%d  %s  %5.1fs\n",
                smp, ncol(sce), median(rho), mean(rho), quantile(rho, .9), length(bg_only), zmode, el))
    rows[[smp]] <- data.frame(sample_id = smp, n_cells = ncol(sce),
                              rho_median = round(median(rho), 5), rho_mean = round(mean(rho), 5),
                              rho_p10 = round(quantile(rho, .10), 5), rho_p90 = round(quantile(rho, .90), 5),
                              n_background_droplets = length(bg_only), z_mode = zmode, seconds = el)
    percell[[smp]] <- data.frame(cell = paste0(smp, "_", colnames(sce)), sample_id = smp,
                                 decontX_contamination = as.numeric(rho))
  }
}
s <- do.call(rbind, rows); p <- do.call(rbind, percell)
stopifnot(!any(duplicated(p$cell)))
write.csv(s, file.path(OUT, "15_ambient_rho_summary.csv"), row.names = FALSE)
write.csv(p, gzfile(file.path(OUT, "15_ambient_rho_percell.csv.gz")), row.names = FALSE)
cat("\n"); print(s, row.names = FALSE)
cat(sprintf("\n跨样本 rho_median 范围: %.4f - %.4f\n", min(s$rho_median), max(s$rho_median)))
cat("STEP15_PASS\n")
