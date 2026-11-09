#!/usr/bin/env Rscript
# =============================================================================
# P51 — Lijian Wu mouse bulk RNA-seq: control vs Crtc2 line x un/TNF/PGE2/adenosine, n=3
# Step 4: DESeq2 differential expression + PCA
#
# Design: 8 groups x 3 replicates, all 24 libraries on ONE run (flowcell 2533TCLT4,
#   lane 3) -> no batch confounding -> one ~Group model; the 6 client contrasts
#   (registration form Section 3G, verbatim, nothing added) are extracted from it.
#   Same design/choice as this client's previous project (proj19, delivered 0724-0729).
# Direction: log2FC > 0 = higher in the FIRST group named (form: "A vs. B, B is baseline").
# Template: proj30 4_run_DE_PCA.R (newest replicated run) + /rnaseq skill mouse filters.
#
# Run: /Work_bio/gao/configs/.conda/envs/DE_R45/bin/Rscript scripts/4_run_DE_PCA.R
#      (never `conda run -n DE_R45`: two envs share that name -> wrong interpreter)
# =============================================================================
suppressPackageStartupMessages({
  library(DESeq2); library(ashr); library(dplyr); library(readr); library(tibble)
  library(ggplot2); library(ggrepel); library(pheatmap); library(RColorBrewer)
})

setwd(dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]))))

TODAY      <- format(Sys.Date(), "%Y%m%d")
DATA_DIR   <- paste0("../Data_Analysis_", TODAY)
DE_DIR     <- file.path(DATA_DIR, "DE_PCA_Results")
READS_DIR  <- file.path(DATA_DIR, "Reads")
RDS_DIR    <- "rds_cache"          # scripts/rds_cache — NEVER inside DATA_DIR
META_FILE  <- "0_sample_table.csv"
COUNT_FILE <- "../output_results/star_salmon/salmon.merged.gene_counts.tsv"
TPM_FILE   <- "../output_results/star_salmon/salmon.merged.gene_tpm.tsv"
for (d in c(DATA_DIR, DE_DIR, READS_DIR, RDS_DIR))
  dir.create(d, showWarnings = FALSE, recursive = TRUE)

PADJ_CUT <- 0.05
LFC_CUT  <- log2(1.2)              # 0.263 -> >= 1.2-fold

# ---- 0. Design -------------------------------------------------------------
TREATS <- c("un", "TNF", "PGE2", "adenosine")
GENOS  <- c("control", "Crtc2")
GROUP_LEVELS <- as.vector(t(outer(GENOS, TREATS, paste, sep = "-")))  # control-un first
# 6 contrasts exactly as in form Section 3G (numerator, denominator=baseline)
CONTRASTS <- list(
  `Crtc2-un_vs_control-un`               = c("Crtc2-un",        "control-un"),
  `Crtc2-TNF_vs_control-TNF`             = c("Crtc2-TNF",       "control-TNF"),
  `Crtc2-PGE2_vs_control-PGE2`           = c("Crtc2-PGE2",      "control-PGE2"),
  `Crtc2-adenosine_vs_control-adenosine` = c("Crtc2-adenosine", "control-adenosine"),
  `control-PGE2_vs_control-un`           = c("control-PGE2",    "control-un"),
  `control-adenosine_vs_control-un`      = c("control-adenosine", "control-un")
)
# V21 paired palette (select_palette(8, relationship="paired")): hue = treatment,
# light = control, dark = Crtc2.
PAL <- c("#F3BB98", "#CE641F", "#819DDD", "#3B5DA4", "#E87878", "#BC3535", "#79CDB8", "#08634F")
GROUP_COLORS <- setNames(PAL, c(rbind(paste0("control-", TREATS), paste0("Crtc2-", TREATS))))
stopifnot(setequal(names(GROUP_COLORS), GROUP_LEVELS))

# ---- 1. Metadata -----------------------------------------------------------
meta <- read_csv(META_FILE, show_col_types = FALSE, col_types = cols(.default = "c")) %>%
  transmute(sample_id, client_sample_name, Group = factor(group, levels = GROUP_LEVELS),
            genotype, treatment, replicate, raw_read_pairs = as.numeric(raw_read_pairs))
stopifnot("Unexpected group label" = !any(is.na(meta$Group)))
stopifnot("Need 3 replicates in every group (full run only)" = all(table(meta$Group) == 3))
print(table(meta$Group))

# ---- 2. Counts -------------------------------------------------------------
raw <- read_tsv(COUNT_FILE, show_col_types = FALSE)
gene_info <- raw %>% transmute(gene_id = gene_id, gene_name = gene_name)
stopifnot("Count matrix lacks samples" = all(meta$sample_id %in% colnames(raw)))
counts_mat <- raw %>% select(all_of(meta$sample_id)) %>% as.matrix() %>% round()
rownames(counts_mat) <- gene_info$gene_id
stopifnot(identical(colnames(counts_mat), meta$sample_id))
cat("Genes in salmon matrix:", nrow(counts_mat), "\n")

reads_summary <- tibble(sample_id = meta$sample_id, client_sample_name = meta$client_sample_name,
                        Group = as.character(meta$Group), raw_read_pairs = meta$raw_read_pairs,
                        assigned_reads = colSums(counts_mat),
                        genes_detected = colSums(counts_mat >= 1))
write_csv(reads_summary, file.path(READS_DIR, "library_metrics_per_sample.csv"))

stopifnot(file.exists(COUNT_FILE), file.exists(TPM_FILE))
file.copy(COUNT_FILE, file.path(READS_DIR, "All_sample_gene_counts.tsv"), overwrite = TRUE)
file.copy(TPM_FILE,   file.path(READS_DIR, "All_sample_gene_tpm.tsv"),    overwrite = TRUE)
stopifnot(file.exists(file.path(READS_DIR, "All_sample_gene_counts.tsv")),
          file.exists(file.path(READS_DIR, "All_sample_gene_tpm.tsv")))

# ---- 2a. Step A — ribosomal AND protein-coding-like (mouse regex, /rnaseq skill) ----
ribo_pattern     <- "^Rpl|^Rps|^Mrpl|^Mrps|^Rplp|^Rpsa"
not_ribo_pattern <- "^Rps6k|^Rps19bp1$|^Rps23rg1$"   # look-alikes that are NOT ribosomal
is_ribo <- grepl(ribo_pattern, gene_info$gene_name, ignore.case = TRUE) &
          !grepl(not_ribo_pattern, gene_info$gene_name, ignore.case = TRUE)
non_ribosomal <- !is_ribo
noncoding_pattern <- paste0("^mt-", "|^Snord|^Snora|^Rnu", "|^Malat1|^Neat1|^Xist",
                            "|^Mir[0-9]|^let-", "|^Linc|^Gm[0-9]", "|pseudogene|antisense")
protein_coding_like <- !grepl(noncoding_pattern, gene_info$gene_name, ignore.case = TRUE)
regex_filter <- non_ribosomal & protein_coding_like       # AND, never if/else (0831 incident)
lookalike <- grepl(not_ribo_pattern, gene_info$gene_name, ignore.case = TRUE)
stopifnot("ribo look-alikes must survive Step A" = all(regex_filter[lookalike & protein_coding_like]))
cat("Step A:", nrow(counts_mat), "->", sum(regex_filter),
    "(ribosomal removed:", sum(is_ribo), "| look-alikes kept:", sum(lookalike), ")\n")
counts_pc <- counts_mat[regex_filter, , drop = FALSE]
gene_pc   <- gene_info[regex_filter, ]

# ---- 2b. Step B — low-expression filter ------------------------------------
n_samples <- nrow(meta)
# Rule v3 (user directive 2026-11-08: "只要比对的两组出现3个样，或者总样品中30%中出现即可"):
#   a gene is testable in a contrast if >=10 counts in >=3 of that contrast's 6 samples
#   (DESeq2 vignette "smallest group size" rule applied to the compared pair; Love's per-contrast
#   advice on support.bioconductor.org/p/109669) OR >=10 counts in >=30% of ALL samples (cf. GTEx
#   20%-of-samples rule). Global pre-filter = union over contrasts; per-contrast restriction +
#   BH happen in section 5. The house n-2 rule (>=22/24) removed condition-specific genes incl.
#   Cxcl1/Cxcl2/Cxcl10/Ccl2/Il6 -> archived as v1 in superseded_for_comparison/DE_enrichment/.
MIN_CNT <- 10L; PAIR_MIN <- 3L; GLOBAL_PROP <- 0.30
grp_of <- as.character(meta$Group[match(colnames(counts_pc), meta$sample_id)])
stopifnot(!anyNA(grp_of), length(grp_of) == n_samples)
global_ok <- rowSums(counts_pc >= MIN_CNT) >= ceiling(GLOBAL_PROP * n_samples)
pair_ok <- sapply(CONTRASTS, function(cc)
  rowSums(counts_pc[, grp_of %in% cc, drop = FALSE] >= MIN_CNT) >= PAIR_MIN)
stopifnot(all(sapply(CONTRASTS, function(cc) sum(grp_of %in% cc)) == 6))
eligible <- pair_ok | global_ok                      # genes x contrasts, used in section 5
keep <- rowSums(eligible) > 0
keep_n2 <- rowSums(counts_pc >= MIN_CNT) >= (n_samples - 2)       # old rule, audit only
stopifnot(all(keep[keep_n2]))
counts_filt <- counts_pc[keep, , drop = FALSE]
gene_filt   <- gene_pc[keep, ]
eligible <- eligible[keep, , drop = FALSE]
cat("Step B (pair>=3 OR >=30% samples):", nrow(counts_pc), "->", nrow(counts_filt),
    "| old n-2 rule would keep", sum(keep_n2), "\n")
# filter-impact audit (§2A): genes kept now but excluded by the old rule, with groups passing
audit <- data.frame(gene_id = rownames(counts_pc)[keep & !keep_n2],
                    gene_name = gene_pc$gene_name[keep & !keep_n2],
                    contrasts_eligible = apply(eligible[!keep_n2[keep], , drop = FALSE], 1,
                                           function(z) paste(colnames(eligible)[z], collapse = ";")))
write.csv(audit, file.path(RDS_DIR, "filter_audit_added_vs_n-2.csv"), row.names = FALSE)
saveRDS(list(n_original = nrow(counts_mat), n_after_regex = sum(regex_filter),
             n_ribo_removed = sum(is_ribo), n_final = nrow(counts_filt),
             n_samples = n_samples, low_count_min = 10L,
             filter_rule = ">=10 counts in >=3 of the 6 samples of the compared pair OR in >=30% of all samples",
             global_prop_min_samples = ceiling(GLOBAL_PROP * n_samples),
             n_eligible_per_contrast = colSums(eligible),
             n_kept_old_n2_rule = sum(keep_n2), n_added_vs_n2 = nrow(audit)),
        file.path(RDS_DIR, "filter_stats.rds"))
id2name <- setNames(gene_filt$gene_name, gene_filt$gene_id)

# ---- 3. DESeq2 -------------------------------------------------------------
dds <- DESeqDataSetFromMatrix(counts_filt, colData = as.data.frame(meta), design = ~ Group)
dds <- DESeq(dds)
vsd <- vst(dds, blind = FALSE)

# ---- 4. PCA ----------------------------------------------------------------
pca_data   <- plotPCA(vsd, intgroup = "Group", returnData = TRUE)
percentVar <- round(100 * attr(pca_data, "percentVar"))
pca_data$name <- rownames(pca_data)
pca_data$genotype  <- meta$genotype[match(pca_data$name, meta$sample_id)]
pca_data$treatment <- meta$treatment[match(pca_data$name, meta$sample_id)]
dist_mat   <- as.matrix(dist(pca_data[, c("PC1", "PC2")]))
same_group <- outer(pca_data$Group, pca_data$Group, "=="); diag(same_group) <- FALSE
within_dist  <- mean(dist_mat[same_group])
between_dist <- mean(dist_mat[!same_group & row(dist_mat) != col(dist_mat)])
separation_ratio <- between_dist / within_dist
cent <- pca_data %>% group_by(Group) %>% summarise(PC1 = mean(PC1), PC2 = mean(PC2), .groups = "drop")
cent_d <- as.matrix(dist(as.matrix(cent[, c("PC1", "PC2")])))
dimnames(cent_d) <- list(as.character(cent$Group), as.character(cent$Group))
dnn <- dist_mat; diag(dnn) <- Inf
nn  <- colnames(dnn)[apply(dnn, 1, which.min)]
nn_tbl <- tibble(sample_id = rownames(dnn), nearest = nn,
                 same_group = as.character(pca_data$Group[match(rownames(dnn), pca_data$name)]) ==
                              as.character(pca_data$Group[match(nn, pca_data$name)]))
cat("\nNearest-neighbour check (FALSE = suspect):\n"); print(as.data.frame(nn_tbl))
# PCs 1-4 too: with 2 crossed factors the genotype axis may sit on PC3+
pca_full <- prcomp(t(assay(vsd)[head(order(-matrixStats::rowVars(assay(vsd))), 500), ]))
pv4 <- round(100 * (pca_full$sdev^2 / sum(pca_full$sdev^2))[1:4], 1)
saveRDS(list(percentVar = percentVar, percentVar4 = pv4, within_dist = within_dist,
             between_dist = between_dist, separation_ratio = separation_ratio,
             pca_data = pca_data, centroid_dist = cent_d, nn_tbl = nn_tbl,
             pc_scores = pca_full$x[, 1:4]),
        file.path(RDS_DIR, "pca_summary.rds"))

p <- ggplot(pca_data, aes(PC1, PC2, color = Group, shape = genotype)) +
  geom_point(size = 4) +
  geom_text_repel(aes(label = name), color = "black", size = 2.6, max.overlaps = Inf, seed = 42) +
  scale_color_manual(values = GROUP_COLORS, breaks = GROUP_LEVELS) +
  scale_shape_manual(values = c(control = 16, Crtc2 = 17)) +
  labs(x = sprintf("PC1: %d%% variance", percentVar[1]),
       y = sprintf("PC2: %d%% variance", percentVar[2]),
       title = sprintf("PCA - all %d samples (VST, top 500 variable genes)", n_samples)) +
  theme_bw(base_size = 12) + theme(legend.position = "right")
ggsave(file.path(DE_DIR, "PCA_AllSamples_Overview.pdf"), p, width = 9.5, height = 6.5)
stopifnot("PCA is missing samples!" = nrow(pca_data) == n_samples)
write_csv(as_tibble(pca_data), file.path(DE_DIR, "PCA_coordinates_overview.csv"))

sampleDists <- dist(t(assay(vsd))); sdm <- as.matrix(sampleDists)
ann_all <- data.frame(Group = as.character(meta$Group[match(colnames(sdm), meta$sample_id)]),
                      row.names = colnames(sdm))
stopifnot("Sample-distance annotation mismatch!" =
  all(ann_all$Group == unname(setNames(as.character(meta$Group), meta$sample_id)[colnames(sdm)])))
pdf(file.path(DE_DIR, "SampleDistance_Heatmap.pdf"), width = 10, height = 9)
pheatmap(sdm, clustering_distance_rows = sampleDists, clustering_distance_cols = sampleDists,
         annotation_col = ann_all, annotation_colors = list(Group = GROUP_COLORS),
         col = colorRampPalette(rev(brewer.pal(9, "Blues")))(255),
         main = "Sample-to-sample distance (VST)")
dev.off()

# ---- 4b. Crtc2 expression itself (the genotype-defining gene; proj19 lesson:
#          always look at whether the target gene behaves as the genotype implies) ----
crtc2_id <- gene_info$gene_id[gene_info$gene_name == "Crtc2"]
nc_all <- counts(estimateSizeFactors(DESeqDataSetFromMatrix(counts_mat, as.data.frame(meta), ~Group)),
                 normalized = TRUE)
crtc2_tbl <- tibble(sample_id = meta$sample_id, Group = as.character(meta$Group),
                    Crtc2_norm_count = as.numeric(nc_all[crtc2_id[1], meta$sample_id]))
write_csv(crtc2_tbl, file.path(RDS_DIR, "crtc2_expression_internal.csv"))
saveRDS(crtc2_tbl, file.path(RDS_DIR, "crtc2_expression.rds"))
print(crtc2_tbl %>% group_by(Group) %>% summarise(mean = mean(Crtc2_norm_count)))

# ---- 5. Contrasts ----------------------------------------------------------
res_list <- list(); deg_counts <- list()
for (cname in names(CONTRASTS)) {
  num <- CONTRASTS[[cname]][1]; den <- CONTRASTS[[cname]][2]
  cat("\n=== Contrast:", cname, "===\n")
  # explicit per-contrast filter => independentFiltering = FALSE, BH within eligible genes
  res <- results(dds, contrast = c("Group", num, den), alpha = PADJ_CUT, independentFiltering = FALSE)
  res <- lfcShrink(dds, contrast = c("Group", num, den), res = res, type = "ashr")
  stopifnot(identical(rownames(res), rownames(eligible)))
  res <- res[eligible[, cname], ]
  res$padj <- p.adjust(res$pvalue, method = "BH")
  smp <- c(meta$sample_id[meta$Group == num], meta$sample_id[meta$Group == den])
  ncnt <- counts(dds, normalized = TRUE)[, smp, drop = FALSE]
  colnames(ncnt) <- paste0("normcount_", smp)
  df <- as.data.frame(res) %>% rownames_to_column("gene_id") %>%
    mutate(gene_name = unname(id2name[gene_id])) %>%
    select(gene_id, gene_name, baseMean, log2FoldChange, lfcSE, pvalue, padj) %>%
    left_join(as.data.frame(ncnt) %>% rownames_to_column("gene_id"), by = "gene_id") %>%
    arrange(padj, desc(abs(log2FoldChange)))
  df$significant <- !is.na(df$padj) & df$padj <= PADJ_CUT & abs(df$log2FoldChange) >= LFC_CUT
  df$direction <- ifelse(!df$significant, "NS", ifelse(df$log2FoldChange > 0, "UP", "DOWN"))

  comp_dir <- file.path(DE_DIR, cname); dir.create(comp_dir, showWarnings = FALSE, recursive = TRUE)
  write_csv(df, file.path(comp_dir, paste0(cname, "_ALL_genes.csv")))
  write_csv(filter(df, significant), file.path(comp_dir, paste0(cname, "_DEG_significant.csv")))
  n_up <- sum(df$direction == "UP"); n_dn <- sum(df$direction == "DOWN")
  cat(sprintf("  UP=%d  DOWN=%d  TOTAL=%d\n", n_up, n_dn, n_up + n_dn))
  sig1 <- !is.na(df$padj) & df$padj <= PADJ_CUT & abs(df$log2FoldChange) >= 1   # form 3H template default
  deg_counts[[cname]] <- tibble(Comparison = cname, Numerator = num, Baseline = den,
                                Up = n_up, Down = n_dn, Total = n_up + n_dn,
                                Up_lfc1 = sum(sig1 & df$log2FoldChange > 0),
                                Down_lfc1 = sum(sig1 & df$log2FoldChange < 0),
                                Genes_tested = nrow(df))

  # volcano — p underflow to 0 => Inf; cap + mark (proj27 label pile-up)
  df$neglog10p <- -log10(df$pvalue)
  finite_max <- max(df$neglog10p[is.finite(df$neglog10p)], na.rm = TRUE)
  y_cap <- ceiling(finite_max / 10) * 10
  df$capped <- is.infinite(df$neglog10p); df$neglog10p[df$capped] <- y_cap
  n_capped <- sum(df$capped, na.rm = TRUE)
  lab_df <- bind_rows(
    filter(df, significant) %>% arrange(padj, desc(abs(log2FoldChange))) %>% head(6),
    filter(df, significant, log2FoldChange > 0) %>% arrange(desc(log2FoldChange)) %>% head(6),
    filter(df, significant, log2FoldChange < 0) %>% arrange(log2FoldChange) %>% head(6)
  ) %>% distinct(gene_id, .keep_all = TRUE)
  sub <- sprintf("padj <= %.2f & |log2FC| >= %.3f   (UP %d / DOWN %d); log2FC > 0 = higher in %s",
                 PADJ_CUT, LFC_CUT, n_up, n_dn, num)
  if (n_capped > 0)
    sub <- paste0(sub, sprintf("\n%d gene(s) with p below double precision shown as triangles at the y-axis cap (%d)",
                               n_capped, y_cap))
  v <- ggplot(df, aes(log2FoldChange, neglog10p, color = direction)) +
    geom_point(aes(shape = capped), alpha = 0.6, size = 1) +
    scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 17), guide = "none") +
    scale_color_manual(values = c(UP = "#BC3535", DOWN = "#3B5DA4", NS = "grey70")) +
    geom_vline(xintercept = c(-LFC_CUT, LFC_CUT), linetype = "dashed", color = "grey40") +
    geom_text_repel(data = lab_df, aes(label = gene_name), color = "black", size = 2.7,
                    max.overlaps = Inf, box.padding = 0.55, point.padding = 0.2,
                    min.segment.length = 0, segment.color = "grey50", segment.size = 0.3, seed = 42) +
    coord_cartesian(ylim = c(0, y_cap * 1.22)) +
    labs(title = paste(num, "vs", den), subtitle = sub, x = "log2 fold change", y = "-log10(p-value)") +
    theme_bw(base_size = 12) + theme(legend.position = "bottom", plot.subtitle = element_text(size = 8))
  stopifnot("Volcano has non-finite y values!" = all(is.finite(df$neglog10p[!is.na(df$pvalue)])))
  ggsave(file.path(comp_dir, paste0(cname, "_volcano.pdf")), v, width = 7.5, height = 7)

  top50 <- head(filter(df, significant)$gene_id, 50)
  if (length(top50) >= 2) {
    smp_comp <- c(meta$sample_id[meta$Group == num], meta$sample_id[meta$Group == den])
    mat <- assay(vsd)[top50, smp_comp, drop = FALSE]
    rn <- ifelse(is.na(id2name[top50]) | id2name[top50] == "", top50, id2name[top50])
    rownames(mat) <- make.unique(rn); mat <- t(scale(t(mat)))
    ann_df <- data.frame(Group = as.character(meta$Group[match(smp_comp, meta$sample_id)]),
                         row.names = smp_comp)
    expected_group <- setNames(as.character(meta$Group), meta$sample_id)[smp_comp]
    stopifnot("Heatmap annotation Group mismatch!" =
                all(as.character(ann_df$Group) == unname(expected_group)))
    pdf(file.path(comp_dir, paste0(cname, "_top50_heatmap.pdf")), width = 7, height = 9)
    pheatmap(mat, annotation_col = ann_df, annotation_colors = list(Group = GROUP_COLORS[c(num, den)]),
             show_rownames = TRUE, fontsize_row = 6, cluster_cols = TRUE,
             main = paste("Top 50 DEGs -", num, "vs", den))
    dev.off()
  }
  res_list[[cname]] <- df
}
deg_summary <- bind_rows(deg_counts)
write_csv(deg_summary, file.path(DE_DIR, "DEG_summary_all_comparisons.csv"))
cat("\n=== DEG summary ===\n"); print(as.data.frame(deg_summary))

norm_counts <- counts(dds, normalized = TRUE) %>% as.data.frame() %>% rownames_to_column("gene_id") %>%
  mutate(gene_name = unname(id2name[gene_id])) %>% relocate(gene_name, .after = gene_id)
write_csv(norm_counts, file.path(DE_DIR, "normalized_counts_all_samples.csv"))

# ---- 6. Cache objects for scripts 5/6/8 (NOT in the client folder) ---------
saveRDS(res_list, file.path(RDS_DIR, "res_list.rds")); saveRDS(vsd, file.path(RDS_DIR, "vsd.rds"))
saveRDS(dds, file.path(RDS_DIR, "dds.rds")); saveRDS(meta, file.path(RDS_DIR, "meta.rds"))
saveRDS(counts_mat, file.path(RDS_DIR, "counts_raw.rds")); saveRDS(deg_summary, file.path(RDS_DIR, "deg_summary.rds"))
saveRDS(CONTRASTS, file.path(RDS_DIR, "contrasts.rds")); saveRDS(reads_summary, file.path(RDS_DIR, "reads_summary.rds"))
saveRDS(GROUP_COLORS, file.path(RDS_DIR, "group_colors.rds"))

# ---- 7. Annotation + MultiQC into the client folder ------------------------
annot_files <- sort(list.files("/Work_bio/references/Genes",
  pattern = "^mouse_Gene_annotation_.*\\.xlsx$", full.names = TRUE), decreasing = TRUE)
if (length(annot_files) > 0)
  file.copy(annot_files[1], file.path(DATA_DIR, basename(annot_files[1])), overwrite = TRUE)
QC_DEST <- file.path(DATA_DIR, "QC"); dir.create(QC_DEST, showWarnings = FALSE, recursive = TRUE)
if (dir.exists("../output_results/multiqc"))
  file.copy("../output_results/multiqc", QC_DEST, recursive = TRUE, overwrite = TRUE)
cat("\nDone. Output:", normalizePath(DATA_DIR), "\n")
