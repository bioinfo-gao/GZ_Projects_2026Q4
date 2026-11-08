#!/usr/bin/env Rscript
# =============================================================================
# P51 — automated pre-delivery audit (adapted from proj30 8_pre_delivery_audit.R).
#   A. IS WHAT IS THERE CORRECT — report numbers vs delivered data, plus INDEPENDENT re-derivations
#      (low-count eligibility from the delivered count matrix; BH padj recomputed per table).
#   B. IS EVERYTHING THERE      — FIXED manifest, not paths scraped from the report.
# It cannot replace looking at the figures. Exit 0 = all passed; any FAIL exits 1.
# Run: /Work_bio/gao/configs/.conda/envs/DE_R45/bin/Rscript 8_pre_delivery_audit.R [<delivery_dir>]
#      (no argument = audit Data_Analysis_*; with a delivery dir = audit that folder instead)
# =============================================================================
suppressPackageStartupMessages({ library(dplyr); library(readr); library(stringr); library(tibble) })
setwd(dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]))))
args <- commandArgs(trailingOnly = TRUE)
data_dirs <- sort(list.dirs("..", full.names = TRUE, recursive = FALSE), decreasing = TRUE)
data_dirs <- data_dirs[grepl("/Data_Analysis_[0-9]{8}$", data_dirs)]
ROOT <- if (length(args)) args[1] else data_dirs[1]
stopifnot(dir.exists(ROOT))
DE_DIR <- file.path(ROOT, "DE_PCA_Results"); ENR_DIR <- file.path(ROOT, "Enrichment"); RDS_DIR <- "rds_cache"
cat("Auditing:", normalizePath(ROOT), "\n")

meta <- readRDS(file.path(RDS_DIR, "meta.rds")); deg <- readRDS(file.path(RDS_DIR, "deg_summary.rds"))
pca <- readRDS(file.path(RDS_DIR, "pca_summary.rds")); reads <- readRDS(file.path(RDS_DIR, "reads_summary.rds"))
fs <- readRDS(file.path(RDS_DIR, "filter_stats.rds")); CON <- readRDS(file.path(RDS_DIR, "contrasts.rds"))

REPORT <- list.files(ROOT, "\\.md$", full.names = TRUE)
fails <- c(); passes <- 0L
check <- function(ok, label, detail = "") {
  if (isTRUE(ok)) { passes <<- passes + 1L; cat(sprintf("  PASS  %s\n", label)) }
  else { fails <<- c(fails, paste0(label, if (nzchar(detail)) paste0(" -- ", detail) else ""))
         cat(sprintf("  FAIL  %s %s\n", label, detail)) } }
check(length(REPORT) == 1, "exactly one report (.md) in the audited folder", paste(basename(REPORT), collapse = ", "))
rep_txt <- readLines(REPORT[1], warn = FALSE); rep_one <- paste(rep_txt, collapse = "\n")
fmt <- function(x) formatC(as.numeric(x), format = "d", big.mark = ",")
has <- function(x) grepl(x, rep_one, fixed = TRUE)

cat("=== 1. sample completeness ===\n")
ov <- read_csv(file.path(DE_DIR, "PCA_coordinates_overview.csv"), show_col_types = FALSE)
check(nrow(ov) == 24 && setequal(ov$name, meta$sample_id), "PCA coordinates contain all 24 samples")
check(all(table(meta$Group) == 3) && length(unique(meta$Group)) == 8, "8 groups x 3 replicates")
st <- read_csv("0_sample_table.csv", show_col_types = FALSE)
check(setequal(st$sample_id, meta$sample_id) && all(st$paired_on_disk == "yes"), "SSOT sample table = meta, all paired")
check(all(meta$client_sample_name == st$client_sample_name[match(meta$sample_id, st$sample_id)]), "client names carried unchanged")

cat("\n=== 2. mislabel screen ===\n")
bad <- pca$nn_tbl %>% filter(!same_group)
if (nrow(bad) == 0) check(TRUE, "every nearest neighbour in own group") else {
  cat(sprintf("  NOTE  %s\n", paste(bad$sample_id, "->", bad$nearest, collapse = "; ")))
  check(all(sapply(bad$sample_id, has)) && has("nearest neighbour"), "report names the off-group nearest neighbours") }

cat("\n=== 3. DEG tables: summary, independent eligibility, BH recompute ===\n")
cnt <- read_tsv(file.path(ROOT, "Reads", "All_sample_gene_counts.tsv"), show_col_types = FALSE)
cm <- round(as.matrix(cnt[, meta$sample_id])); rownames(cm) <- cnt$gene_id   # DESeq2 + filter use rounded estimated counts
grp <- as.character(meta$Group)
glob_ok <- rowSums(cm >= 10) >= ceiling(0.3 * 24)
ribo <- "^(Rpl|Rps|Mrpl|Mrps|Rplp|Rpsa)"; notribo <- "^(Rps6k|Rps19bp1$|Rps23rg1$)"
for (i in seq_len(nrow(deg))) { r <- deg[i, ]; k <- r$Comparison
  dd <- read_csv(file.path(DE_DIR, k, paste0(k, "_DEG_significant.csv")), show_col_types = FALSE)
  al <- read_csv(file.path(DE_DIR, k, paste0(k, "_ALL_genes.csv")), show_col_types = FALSE)
  check(sum(dd$log2FoldChange > 0) == r$Up && sum(dd$log2FoldChange < 0) == r$Down && nrow(dd) == r$Total,
        paste("DEG table matches summary:", k))
  check(nrow(al) == r$Genes_tested, paste("genes tested matches:", k), sprintf("%d vs %d", nrow(al), r$Genes_tested))
  pair_ok <- rowSums(cm[, grp %in% CON[[k]], drop = FALSE] >= 10) >= 3
  elig <- rownames(cm)[(pair_ok | glob_ok)]
  check(all(al$gene_id %in% elig), paste("every tested gene meets the stated low-count rule (independent):", k),
        sprintf("%d violate", sum(!al$gene_id %in% elig)))
  check(isTRUE(all.equal(p.adjust(al$pvalue, "BH"), al$padj, tolerance = 1e-8)), paste("padj = BH within the tested set:", k))
  sig <- !is.na(al$padj) & al$padj <= 0.05 & abs(al$log2FoldChange) >= log2(1.2)
  check(sum(sig) == r$Total, paste("significance rule re-derived:", k))
  check(!any(grepl(ribo, dd$gene_name) & !grepl(notribo, dd$gene_name)), paste("no ribosomal genes among DEGs:", k))
  check("stat" %in% names(al) && all(is.finite(al$stat)), paste("Wald stat column present (GSEA rank):", k)) }
check(fs$n_final == length(unique(unlist(lapply(names(CON), function(k)
  read_csv(file.path(DE_DIR, k, paste0(k, "_ALL_genes.csv")), show_col_types = FALSE)$gene_id)))),
  "global gene count = union of per-comparison tested sets")

cat("\n=== 4. report text vs data ===\n")
for (i in seq_len(nrow(deg))) { r <- deg[i, ]; tag <- gsub("_vs_", " vs ", r$Comparison)
  rows <- grep(paste0("^\\| ", tag, " \\|"), rep_txt, value = TRUE)
  ok <- any(sapply(rows, function(row) { v <- suppressWarnings(as.numeric(gsub("[^0-9.]", "", trimws(str_split(row, "\\|")[[1]]))))
    all(c(r$Up, r$Down, r$Total) %in% v) }))
  check(ok, paste("report DEG row matches data:", tag)) }
check(all(sapply(meta$sample_id, function(s) any(grepl(paste0("^\\| ", s, " \\|"), rep_txt)))), "report lists every sample")
check(all(sapply(reads$raw_read_pairs, function(x) has(fmt(x)))), "report quotes read pairs for every sample")
check(has(sprintf("%.1f", pca$separation_ratio)), "report quotes PCA separation ratio")
check(has(fmt(fs$n_final)), "report quotes number of genes kept")
al1 <- read_csv(file.path(DE_DIR, names(CON)[1], paste0(names(CON)[1], "_ALL_genes.csv")), show_col_types = FALSE)
check(has(sprintf("%+.2f", al1$log2FoldChange[al1$gene_name == "Crtc2"])), "report Crtc2 log2FC (comparison 1) = delivered table")
ret <- read_csv(file.path(DE_DIR, "Crtc2_line_treatment_response_retention.csv"), show_col_types = FALSE)
for (j in seq_len(nrow(ret))) check(has(sprintf("%.0f%%", 100 * ret$Slope_Crtc2_vs_control_response[j])),
  paste("report response fraction = delivered CSV:", ret$Treatment[j]))

cat("\n=== 5. metadata discipline (only what the client stated) ===\n")
check(!grepl("knock-?out|knock-?down|overexpress|\\bKO\\b|\\bKD\\b", rep_one, ignore.case = TRUE),
      "report never states how Crtc2 was modified")
check(has("cell line name not specified by client"), "header states the cell line name was not given")
check(has("Mus musculus") && has("Tissue / Cell:"), "header carries species and tissue/cell")
check(!grepl("RAW ?264|NIH ?3T3|MEF|macrophage|fibroblast|T cell line|Jurkat|EL4", rep_one, ignore.case = TRUE),
      "report names no specific cell line / cell type the client did not state")
check(!grepl("long-read|nanopore|pacbio", rep_one, ignore.case = TRUE), "no long-read recommendation")
check(has("Zhen Gao, PhD, Principal Bioinformatics Scientist, Athenomics"), "signature present")

cat("\n=== 6. structure and cross-references ===\n")
for (s in c("## 1. Objectives", "## 2. Core Conclusions", "## 3. Key Findings", "## 4. Sample Information",
            "## 5. Analysis Rationale and Decision Criteria", "## 6. Methods", "## 7. Results",
            "## 8. Conclusions, Limitations and Next Steps", "## 9. Deliverable Files"))
  check(any(rep_txt == s), paste("section present:", s))
check(any(grepl("^### 8.2 Practical reading", rep_txt)) && any(grepl("^### 7.1 Principal Component Analysis", rep_txt)),
      "8.2 and 7.1 present")
refs <- unique(unlist(str_extract_all(rep_one, "[Ss]ection [0-9]+(\\.[0-9]+)?")))
heads <- c(str_match(rep_txt, "^## ([0-9]+)\\.")[, 2], str_match(rep_txt, "^### ([0-9]+\\.[0-9]+)")[, 2]); heads <- heads[!is.na(heads)]
badref <- refs[!sub("^[Ss]ection ", "", refs) %in% heads]
check(length(badref) == 0, "every Section N.M cross-reference resolves", paste(badref, collapse = ", "))
tbl <- rep_txt[grepl("^\\| *:-", rep_txt)]
check(all(grepl("^\\| :--- \\|( :---: \\|)+$", tbl)), "table alignment rows: first :---, rest :---:")
hdr_idx <- which(grepl("^\\| *:-", rep_txt)) - 1
ncol_ok <- sapply(hdr_idx, function(h) { n <- str_count(rep_txt[h], "\\|"); j <- h + 1
  while (j <= length(rep_txt) && grepl("^\\|", rep_txt[j])) { if (str_count(rep_txt[j], "\\|") != n) return(FALSE); j <- j + 1 }; TRUE })
check(all(ncol_ok), "every table row has the same number of cells as its header")

cat("\n=== 7. fixed manifest ===\n")
MAN <- c("Reads/All_sample_gene_counts.tsv", "Reads/All_sample_gene_tpm.tsv", "Reads/library_metrics_per_sample.csv",
         "QC/multiqc", "DE_PCA_Results/PCA_AllSamples_Overview.pdf", "DE_PCA_Results/PCA_coordinates_overview.csv",
         "DE_PCA_Results/SampleDistance_Heatmap.pdf", "DE_PCA_Results/normalized_counts_all_samples.csv",
         "DE_PCA_Results/DEG_summary_all_comparisons.csv", "DE_PCA_Results/Crtc2_line_treatment_response_retention.csv",
         "DE_PCA_Results/CREB_CRTC_target_and_inflammatory_marker_panel.csv",
         "Enrichment/Enrichment_summary_all_comparisons.csv", "mouse_Gene_annotation_20260628.xlsx")
for (p in MAN) check(file.exists(file.path(ROOT, p)), paste("manifest:", p))
for (k in names(CON)) for (sub in c(file.path("DE_PCA_Results", k), file.path("Enrichment", k, c("GO", "KEGG", "GSEA"))))
  check(dir.exists(file.path(ROOT, sub)), paste("manifest:", sub))
sec9 <- rep_txt[seq(which(rep_txt == "## 9. Deliverable Files"), length(rep_txt))]
for (p in MAN) check(any(grepl(p, sec9, fixed = TRUE)), paste("Section 9 lists:", p))
paths <- unique(gsub("`", "", unlist(str_extract_all(rep_one, "`[^`]+`"))))
paths <- paths[grepl("/|\\.(csv|pdf|tsv|xlsx)$", paths) & !grepl("[<>*{]|~ |\\$", paths)]
miss <- paths[!file.exists(file.path(ROOT, paths)) & !file.exists(file.path(DE_DIR, paths))]
check(length(miss) == 0, "every concrete path named in the report exists", paste(miss, collapse = ", "))

cat("\n=== 8. hygiene and figures ===\n")
junk <- list.files(ROOT, "\\.(rds|RData|log)$", recursive = TRUE)
check(length(junk) == 0, "no .rds/.RData/.log", paste(head(junk), collapse = ", "))
png <- list.files(ROOT, "\\.png$", recursive = TRUE); png <- png[!grepl("^QC/", png)]
check(length(png) == 0, "no PNG among our own figures", paste(head(png), collapse = ", "))
pdfs <- c(file.path(DE_DIR, c("PCA_AllSamples_Overview.pdf", "SampleDistance_Heatmap.pdf")),
          unlist(lapply(names(CON), function(k) file.path(DE_DIR, k, paste0(k, c("_volcano.pdf", "_top50_heatmap.pdf"))))))
for (p in pdfs) check(file.exists(p) && file.size(p) > 1000, paste("figure present:", basename(p)))

cat("\n==================== AUDIT SUMMARY ====================\n")
cat(sprintf("passed: %d   failed: %d\n", passes, length(fails)))
if (length(fails)) { cat("FAILURES:\n"); for (f in fails) cat("  - ", f, "\n"); quit(status = 1) }
cat("ALL AUTOMATED CHECKS PASSED. Still required: open the PCA, volcano and heatmap PDFs.\n")
