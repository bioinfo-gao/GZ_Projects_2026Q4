#!/usr/bin/env Rscript
# =============================================================================
# P51 — Lijian Wu mouse RNA-seq. Step 5: GO / KEGG ORA (UP/DOWN/ALL) + GSEA
#   (GO:BP, KEGG, MSigDB Hallmark mouse) for the 6 client contrasts.
# Template: proj30 5_run_enrichment.R, adapted to mouse (org.Mm.eg.db / mmu).
# Change vs template: a FAILED call (e.g. KEGG REST unreachable) is recorded as
#   status=ERROR in the summary instead of being silently skipped — a skipped
#   call previously looked identical to "0 significant pathways".
# Run: /Work_bio/gao/configs/.conda/envs/DE_R45/bin/Rscript scripts/5_run_enrichment.R
# =============================================================================
suppressPackageStartupMessages({
  library(clusterProfiler); library(org.Mm.eg.db); library(enrichplot)
  library(msigdbr); library(dplyr); library(readr); library(ggplot2); library(tibble)
})
setwd(dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]))))

data_dirs <- sort(list.dirs("..", full.names = TRUE, recursive = FALSE), decreasing = TRUE)
data_dirs <- data_dirs[grepl("/Data_Analysis_[0-9]{8}$", data_dirs)]
if (length(data_dirs) == 0) stop("No Data_Analysis_YYYYMMDD folder found.")
DATA_DIR <- data_dirs[1]
ENR_DIR  <- file.path(DATA_DIR, "Enrichment"); RDS_DIR <- "rds_cache"
dir.create(ENR_DIR, showWarnings = FALSE, recursive = TRUE)
cat("Data dir:", DATA_DIR, "\n")
res_list <- readRDS(file.path(RDS_DIR, "res_list.rds"))
ORG <- org.Mm.eg.db; KEGG_ORG <- "mmu"

hallmark <- tryCatch(msigdbr(species = "Mus musculus", collection = "H"),
                     error = function(e) msigdbr(species = "Mus musculus", category = "H")) %>%
  dplyr::select(gs_name, gene_symbol)
stopifnot("Hallmark gene sets empty" = length(unique(hallmark$gs_name)) == 50)

save_dot <- function(obj, path, title, showN = 20) {
  if (is.null(obj) || nrow(as.data.frame(obj)) == 0) return(invisible(NULL))
  n <- min(showN, nrow(as.data.frame(obj)))
  p <- try(dotplot(obj, showCategory = n) + ggtitle(title) +
             theme(axis.text.y = element_text(size = 7)), silent = TRUE)
  if (!inherits(p, "try-error"))
    ggsave(path, p, width = 9, height = max(4, 0.28 * n + 2), limitsize = FALSE)
  invisible(NULL)
}
rows <- list()
rec <- function(cname, analysis, status, n, n_input = NA)
  rows[[length(rows) + 1]] <<- tibble(Comparison = cname, Analysis = analysis,
                                       status = status, N_significant = n, n_input = n_input)
run1 <- function(cname, analysis, expr, csv, dotpath, title, n_input = NA) {
  obj <- tryCatch(expr, error = function(e) e)
  if (inherits(obj, "error")) {
    cat("  !! ", analysis, " ERROR: ", conditionMessage(obj), "\n", sep = "")
    rec(cname, analysis, "ERROR", NA_integer_, n_input); return(invisible(NULL))
  }
  tab <- if (is.null(obj)) data.frame() else as.data.frame(obj)
  write_csv(tab, csv); save_dot(obj, dotpath, title)
  rec(cname, analysis, "ok", nrow(tab), n_input)
}

for (cname in names(res_list)) {
  cat("\n########## ", cname, " ##########\n")
  df <- res_list[[cname]]
  go_dir <- file.path(ENR_DIR, cname, "GO"); kegg_dir <- file.path(ENR_DIR, cname, "KEGG")
  gsea_dir <- file.path(ENR_DIR, cname, "GSEA")
  for (d in c(go_dir, kegg_dir, gsea_dir)) dir.create(d, recursive = TRUE, showWarnings = FALSE)

  map <- suppressWarnings(bitr(unique(df$gene_name[!is.na(df$gene_name) & df$gene_name != ""]),
                               fromType = "SYMBOL", toType = "ENTREZID", OrgDb = ORG))
  df2 <- df %>% inner_join(map, by = c("gene_name" = "SYMBOL"))
  universe <- unique(df2$ENTREZID)
  stopifnot("ENTREZ mapping collapsed" = length(universe) > 0.7 * nrow(df))
  gene_sets <- list(
    UP   = df2 %>% filter(significant, log2FoldChange > 0) %>% pull(ENTREZID) %>% unique(),
    DOWN = df2 %>% filter(significant, log2FoldChange < 0) %>% pull(ENTREZID) %>% unique(),
    ALL  = df2 %>% filter(significant) %>% pull(ENTREZID) %>% unique())
  cat(sprintf("  ENTREZ-mapped DEGs: UP=%d DOWN=%d ALL=%d (universe %d)\n",
              length(gene_sets$UP), length(gene_sets$DOWN), length(gene_sets$ALL), length(universe)))

  for (ont in c("BP", "MF", "CC")) for (dl in names(gene_sets)) {
    g <- gene_sets[[dl]]; an <- paste0("GO_", ont, "_", dl)
    if (length(g) < 5) { rec(cname, an, "skipped_lt5_genes", NA_integer_, length(g)); next }
    run1(cname, an, enrichGO(g, OrgDb = ORG, universe = universe, ont = ont, keyType = "ENTREZID",
                             pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.2, readable = TRUE),
         file.path(go_dir, sprintf("GO_%s_%s.csv", ont, dl)),
         file.path(go_dir, sprintf("GO_%s_%s_dotplot.pdf", ont, dl)),
         sprintf("%s - GO:%s (%s)", cname, ont, dl), length(g))
  }
  for (dl in names(gene_sets)) {
    g <- gene_sets[[dl]]; an <- paste0("KEGG_", dl)
    if (length(g) < 5) { rec(cname, an, "skipped_lt5_genes", NA_integer_, length(g)); next }
    run1(cname, an, {
      ek <- enrichKEGG(g, organism = KEGG_ORG, universe = universe, pAdjustMethod = "BH",
                       pvalueCutoff = 0.05, qvalueCutoff = 0.2)
      if (!is.null(ek)) ek <- setReadable(ek, OrgDb = ORG, keyType = "ENTREZID"); ek },
      file.path(kegg_dir, sprintf("KEGG_%s.csv", dl)),
      file.path(kegg_dir, sprintf("KEGG_%s_dotplot.pdf", dl)),
      sprintf("%s - KEGG (%s)", cname, dl), length(g))
  }
  # GSEA rank = DESeq2 Wald statistic (user directive 2026-11-08). The house metric
  # sign(LFC) x -log10(p) is so heavy-tailed here (up to 871 genes with |value| > 50) that the
  # weighted ES is dominated by a few extreme genes: Crtc2-PGE2 vs control-PGE2 gave 0 GO:BP terms
  # (min padj 0.063) vs 479 with the Wald stat, same top terms. v3 signlogp version archived.
  stopifnot("Wald stat missing from res_list (rerun 4_run_DE_PCA.R)" = "stat" %in% names(df2))
  rank_df <- df2 %>% filter(!is.na(stat), is.finite(stat)) %>%
    mutate(stat_rank = stat) %>%
    group_by(ENTREZID) %>% slice_max(abs(stat_rank), n = 1, with_ties = FALSE) %>% ungroup()
  gl <- sort(setNames(rank_df$stat_rank, rank_df$ENTREZID), decreasing = TRUE)
  gs <- rank_df %>% group_by(gene_name) %>% slice_max(abs(stat_rank), n = 1, with_ties = FALSE) %>% ungroup()
  gl_sym <- sort(setNames(gs$stat_rank, gs$gene_name), decreasing = TRUE)

  set.seed(42)
  run1(cname, "GSEA_GO_BP", {
    x <- gseGO(gl, OrgDb = ORG, ont = "BP", keyType = "ENTREZID", pvalueCutoff = 0.05,
               pAdjustMethod = "BH", seed = TRUE, verbose = FALSE)
    if (!is.null(x)) x <- setReadable(x, OrgDb = ORG, keyType = "ENTREZID"); x },
    file.path(gsea_dir, "GSEA_GO_BP.csv"), file.path(gsea_dir, "GSEA_GO_BP_dotplot.pdf"),
    sprintf("%s - GSEA GO:BP", cname), length(gl))
  set.seed(42)
  run1(cname, "GSEA_KEGG", {
    x <- gseKEGG(gl, organism = KEGG_ORG, pvalueCutoff = 0.05, pAdjustMethod = "BH",
                 seed = TRUE, verbose = FALSE)
    if (!is.null(x)) x <- setReadable(x, OrgDb = ORG, keyType = "ENTREZID"); x },
    file.path(gsea_dir, "GSEA_KEGG.csv"), file.path(gsea_dir, "GSEA_KEGG_dotplot.pdf"),
    sprintf("%s - GSEA KEGG", cname), length(gl))
  set.seed(42)
  run1(cname, "GSEA_Hallmark",
       GSEA(gl_sym, TERM2GENE = hallmark, pvalueCutoff = 0.05, pAdjustMethod = "BH",
            seed = TRUE, verbose = FALSE),
       file.path(gsea_dir, "GSEA_MSigDB_Hallmark.csv"),
       file.path(gsea_dir, "GSEA_MSigDB_Hallmark_dotplot.pdf"),
       sprintf("%s - GSEA MSigDB Hallmark", cname), length(gl_sym))
}
summ <- bind_rows(rows)
write_csv(summ, file.path(ENR_DIR, "Enrichment_summary_all_comparisons.csv"))
saveRDS(summ, file.path(RDS_DIR, "enrich_summary.rds"))
n_err <- sum(summ$status == "ERROR")
cat("\nDone. Enrichment written to", normalizePath(ENR_DIR), "| ERROR calls:", n_err, "\n")
if (n_err > 0) { print(as.data.frame(filter(summ, status == "ERROR"))); quit(status = 2) }
