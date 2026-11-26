#!/usr/bin/env Rscript
# P47 round1 step 13a：逐样本 scDblFinder。
# 为什么逐样本而不是合并后跑：/scrna §3.2 —— 合并后跑会把样本间差异当成 doublet 信号。
# 为什么用 R 而不是 scanpy 的 scrublet：skill 判定 scDblFinder benchmark 最强；
#   本 env(p47_r) 的 xgboost 已钉 3.3.1.1 并经 .xgbtrain 直接调用验证可训练(step 11)。
# 不搭桥：直接读 staged 的 10x h5，避免 h5ad<->SCE 转置/元数据丢失那一类静默错误。
suppressPackageStartupMessages({library(scDblFinder); library(SingleCellExperiment); library(Seurat); library(Matrix)})
set.seed(20261125)

STAGE <- "/mnt/ex_8T_SSD/projects_2026Q4/47_Liu_Yuan_scRNA/staged_inputs_1125"
OUT   <- "/home/gao/projects_2026Q4/47_Liu_Yuan_scRNA/output_results/round1_descriptive_1125"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
cat("scDblFinder", as.character(packageVersion("scDblFinder")),
    "xgboost", as.character(packageVersion("xgboost")), "\n")

libs <- c("N1-ABE-N2-ABE-N1-CBE-N2-CBE", "N1-Control-N2-Control-CBE-Control-ABE-Control")
res <- list(); summ <- list()
for (lib in libs) {
  base <- file.path(STAGE, lib, "per_sample_outs")
  for (smp in sort(list.dirs(base, recursive = FALSE, full.names = FALSE))) {
    h5 <- file.path(base, smp, "sample_filtered_feature_bc_matrix.h5")
    m  <- Read10X_h5(h5)
    sce <- SingleCellExperiment(list(counts = m))
    t0 <- Sys.time()
    sce <- scDblFinder(sce, verbose = FALSE)
    el <- round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1)
    cls <- sce$scDblFinder.class; scr <- sce$scDblFinder.score
    # 阳性对照式断言：分数必须真的有分布，不能是常数(=训练被静默吞错后的回退分数特征)
    stopifnot(length(unique(round(scr, 6))) > 100)
    cat(sprintf("  %-12s cells=%6d doublet=%5d (%.2f%%)  score sd=%.4f  %5.1fs\n",
                smp, ncol(sce), sum(cls == "doublet"), 100 * mean(cls == "doublet"), sd(scr), el))
    res[[smp]] <- data.frame(cell = paste0(smp, "_", colnames(sce)), sample_id = smp,
                             scDblFinder_score = scr, scDblFinder_class = as.character(cls))
    summ[[smp]] <- data.frame(sample_id = smp, n_cells = ncol(sce),
                              n_doublet = sum(cls == "doublet"),
                              pct_doublet = round(100 * mean(cls == "doublet"), 3),
                              score_sd = round(sd(scr), 4), seconds = el)
  }
}
d <- do.call(rbind, res); s <- do.call(rbind, summ)
stopifnot(!any(duplicated(d$cell)))
write.csv(d, file.path(OUT, "13a_doublet_percell.csv"), row.names = FALSE)
write.csv(s, file.path(OUT, "13a_doublet_summary.csv"), row.names = FALSE)
cat("\n"); print(s, row.names = FALSE)
cat(sprintf("\ntotal cells=%d  total doublets=%d (%.2f%%)\n", nrow(d), sum(d$scDblFinder_class == "doublet"),
            100 * mean(d$scDblFinder_class == "doublet")))
cat("STEP13A_PASS\n")
