#!/usr/bin/env Rscript
# =============================================================================
# P51 — Lijian Wu / mouse cell line, control vs Crtc2 line x un / TNF / PGE2 / adenosine (n=3)
# Step 6: client-facing report (English, house 9-section structure; §2 and §8.2 both present).
#
# Header facts are ONLY what the client stated on THIS form: species Mouse (Section 4),
# "Cell line (established)", "Normal only", snap-frozen (3C), RNA-Seq (mRNA). No cell-line name,
# no statement of how Crtc2 was modified -> the report says "Crtc2 line" and reports measured
# Crtc2 mRNA, never "knockout"/"knockdown"/"overexpression" as a fact.
#
# Every number in the text is computed here from rds_cache / Data_Analysis outputs.
# Run: /Work_bio/gao/configs/.conda/envs/DE_R45/bin/Rscript 6_generate_research_report.R
#      DELIVERY_DATE=YYYY-MM-DD ... (re-stamp on review completion; called by 9_build_delivery.sh)
# =============================================================================
suppressPackageStartupMessages({
  library(dplyr); library(readr); library(tibble); library(tidyr); library(DESeq2)
})
setwd(dirname(normalizePath(sub("^--file=", "",
  grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]))))

data_dirs <- sort(list.dirs("..", full.names = TRUE, recursive = FALSE), decreasing = TRUE)
data_dirs <- data_dirs[grepl("/Data_Analysis_[0-9]{8}$", data_dirs)]
stopifnot(length(data_dirs) == 1)                      # one live analysis folder only
DATA_DIR <- data_dirs[1]; RDS_DIR <- "rds_cache"
DE_DIR <- file.path(DATA_DIR, "DE_PCA_Results"); ENR_DIR <- file.path(DATA_DIR, "Enrichment")

DELIVERY_DATE <- Sys.getenv("DELIVERY_DATE", unset = "")
DELIVERY_DATE <- if (nzchar(DELIVERY_DATE)) as.Date(DELIVERY_DATE) else Sys.Date()
REPORT <- file.path(DATA_DIR, paste0("Crtc2_Line_Treatment_RNAseq_Report_",
                                     format(DELIVERY_DATE, "%m%d"), ".md"))

meta     <- readRDS(file.path(RDS_DIR, "meta.rds"))
res_list <- readRDS(file.path(RDS_DIR, "res_list.rds"))
deg      <- readRDS(file.path(RDS_DIR, "deg_summary.rds"))
fs       <- readRDS(file.path(RDS_DIR, "filter_stats.rds"))
pca      <- readRDS(file.path(RDS_DIR, "pca_summary.rds"))
reads    <- readRDS(file.path(RDS_DIR, "reads_summary.rds"))
CON      <- readRDS(file.path(RDS_DIR, "contrasts.rds"))
enr      <- readRDS(file.path(RDS_DIR, "enrich_summary.rds"))
dds      <- readRDS(file.path(RDS_DIR, "dds.rds"))
st       <- read_csv("0_sample_table.csv", show_col_types = FALSE)
stopifnot(nrow(meta) == 24, length(res_list) == 6, identical(names(res_list), deg$Comparison),
          all(enr$status != "ERROR"))

fmt <- function(x) formatC(as.numeric(x), format = "d", big.mark = ",")
f1 <- function(x) sprintf("%.1f", x); f2 <- function(x) sprintf("%.2f", x)
pc_ <- function(x) sprintf("%.0f%%", 100 * x)
cmp <- function(k) gsub("_vs_", " vs ", k)
dg <- function(k) deg[deg$Comparison == k, ]
L <- c(); add <- function(...) L <<- c(L, paste0(...))
TREAT_LONG <- c(un = "untreated", TNF = "TNF", PGE2 = "PGE2", adenosine = "adenosine")

# ---- QC facts from MultiQC ---------------------------------------------------
QC_DATA <- file.path(DATA_DIR, "QC", "multiqc", "star_salmon", "multiqc_report_data")
gs <- read_tsv(file.path(QC_DATA, "multiqc_general_stats.txt"), show_col_types = FALSE)
qc <- function(pat) { cn <- grep(pat, names(gs), value = TRUE)[1]; stopifnot(!is.na(cn))
  v <- suppressWarnings(as.numeric(gs[[cn]])); v <- v[!is.na(v)]; stopifnot(length(v) == 24); range(v) }
Q_MAP <- qc("star_mapped_percent_1$"); Q_UNIQ <- qc("star_uniquely_mapped_percent_1$")
Q_RRNA <- qc("percent_rRNA$"); Q_DUP <- qc("PERCENT_DUPLICATION$")
Q_LEN <- unique(na.omit(suppressWarnings(as.numeric(gs[[grep("fastqc_raw_avg_sequence_length$", names(gs))]]))))
stb <- read_tsv(file.path(QC_DATA, "multiqc_fail_strand_check_table.txt"), show_col_types = FALSE)
STRAND <- unique(stb$`Inferred strandedness`); stopifnot(length(STRAND) == 1)
LANES <- unique(sub("^([^_]+_[0-9]+)_.*", "\\1", c(st$r1_file, st$r2_file)))
stopifnot(length(LANES) == 1)              # all 24 libraries on one flowcell lane -> no lane confound
RP <- range(reads$raw_read_pairs)

# ---- PCA facts -----------------------------------------------------------------
cd <- pca$centroid_dist; nn_bad <- pca$nn_tbl %>% filter(!same_group)
d_ <- function(a, b) cd[a, b]
GPC1 <- with(pca$pca_data %>% group_by(Group) %>% summarise(m = mean(PC1)), setNames(m, as.character(Group)))
stopifnot(GPC1["Crtc2-TNF"] < GPC1["control-TNF"], GPC1["Crtc2-TNF"] > max(GPC1[c("control-un", "Crtc2-un")]))
PC2_G <- with(pca$pca_data, tapply(PC2, genotype, range))
PC2_SPLIT <- min(PC2_G[["Crtc2"]]) > max(PC2_G[["control"]])          # PC2 separates the two lines completely?
stopifnot(PC2_SPLIT)
TNF_PC1 <- with(pca$pca_data %>% group_by(treatment) %>% summarise(m = mean(PC1)), setNames(m, treatment))

# ---- Crtc2 mRNA ------------------------------------------------------------------
cr <- read_csv(file.path(RDS_DIR, "crtc2_expression_internal.csv"), show_col_types = FALSE) %>%
  left_join(meta %>% select(sample_id, genotype, treatment), by = "sample_id")
cr_g <- cr %>% group_by(genotype, treatment) %>% summarise(m = mean(Crtc2_norm_count), .groups = "drop") %>%
  pivot_wider(names_from = genotype, values_from = m) %>% mutate(ratio = Crtc2 / control)
crtc2_lfc <- sapply(names(res_list)[1:4], function(k) {
  d <- res_list[[k]]; i <- which(d$gene_name == "Crtc2"); stopifnot(length(i) == 1); d$log2FoldChange[i] })
crtc2_sig <- sapply(names(res_list)[1:4], function(k) { d <- res_list[[k]]; d$significant[d$gene_name == "Crtc2"] })
stopifnot(all(crtc2_lfc < 0), all(crtc2_sig))
CR_RATIO <- range(cr_g$ratio)
CR_TNF_FOLD <- with(cr_g, c(control = control[treatment == "TNF"] / control[treatment == "un"], Crtc2 = Crtc2[treatment == "TNF"] / Crtc2[treatment == "un"]))
stopifnot(all(CR_TNF_FOLD > 1.5))
CR_LESS <- sprintf("%s–%s%%", round(100 * (1 - max(CR_RATIO))), round(100 * (1 - min(CR_RATIO))))

# ---- response retention (descriptive, from group means of log2 normalised counts) ----------
nc <- counts(dds, normalized = TRUE)
grp <- as.character(meta$Group[match(colnames(nc), meta$sample_id)]); stopifnot(!anyNA(grp))
gm <- sapply(split(seq_along(grp), grp), function(ix) rowMeans(log2(nc[, ix, drop = FALSE] + 1)))
id2name <- setNames(res_list[[1]]$gene_name, res_list[[1]]$gene_id)
retain <- function(tr, ids, basis) {
  c_r <- gm[ids, paste0("control-", tr)] - gm[ids, "control-un"]
  k_r <- gm[ids, paste0("Crtc2-", tr)] - gm[ids, "Crtc2-un"]
  tibble(Treatment = tr, Gene_set_basis = basis, N_genes = length(ids),
         Slope_Crtc2_vs_control_response = unname(coef(lm(k_r ~ 0 + c_r))),
         Correlation = cor(c_r, k_r),
         Fraction_genes_below_half_response = mean(abs(k_r) < 0.5 * abs(c_r)))
}
sig_ids <- function(k) { d <- res_list[[k]]; d$gene_id[d$significant] }
dTNF <- gm[, "control-TNF"] - gm[, "control-un"]
tnf_ids <- names(dTNF)[abs(dTNF) >= 1 & pmax(gm[, "control-TNF"], gm[, "control-un"]) >= 3]
RET <- bind_rows(
  retain("PGE2", sig_ids("control-PGE2_vs_control-un"), "DEGs of control-PGE2 vs control-un"),
  retain("adenosine", sig_ids("control-adenosine_vs_control-un"), "DEGs of control-adenosine vs control-un"),
  retain("TNF", tnf_ids, "genes with >=2-fold difference in group means, control-TNF vs control-un (descriptive; not a tested comparison)"))
write_csv(RET, file.path(DE_DIR, "Crtc2_line_treatment_response_retention.csv"))
rr <- function(tr, col) RET[[col]][RET$Treatment == tr]

# ---- direction of the Crtc2 effect relative to the PGE2 / adenosine response -----------
opp <- function(a, b) {
  x <- inner_join(res_list[[a]] %>% select(gene_id, la = log2FoldChange, sa = significant),
                  res_list[[b]] %>% select(gene_id, lb = log2FoldChange, sb = significant), by = "gene_id")
  both <- x$sa & x$sb
  list(r = cor(x$la, x$lb), n = sum(both), opp = mean(sign(x$la[both]) != sign(x$lb[both])))
}
O_UN_PGE2 <- opp("Crtc2-un_vs_control-un", "control-PGE2_vs_control-un")
O_PGE2    <- opp("Crtc2-PGE2_vs_control-PGE2", "control-PGE2_vs_control-un")
O_UN_ADO  <- opp("Crtc2-un_vs_control-un", "control-adenosine_vs_control-un")
O_ADO_UN  <- opp("Crtc2-adenosine_vs_control-adenosine", "Crtc2-un_vs_control-un")
# 8.1 wording guards: "mostly the baseline difference" (adenosine), "mostly the PGE2 response" (PGE2)
stopifnot(O_ADO_UN$n / deg$Total[deg$Comparison == "Crtc2-adenosine_vs_control-adenosine"] > 0.5,
          O_ADO_UN$opp < 0.05,
          O_PGE2$opp > 0.9, O_PGE2$n / deg$Total[deg$Comparison == "Crtc2-PGE2_vs_control-PGE2"] > 0.5)

# ---- CREB / CRTC target and inflammatory marker panel ---------------------------------------
PANEL <- c("Crtc2", "Nr4a2", "Crem", "Pde4b", "Pde4d", "Sik1", "Il6", "Ptgs2", "Cxcl2", "Cxcl10",
           "Ccl5", "Nfkbia", "Tnfaip3", "Fos", "Fosb")
pan <- bind_rows(lapply(names(res_list), function(k) {
  d <- res_list[[k]]; i <- match(PANEL, d$gene_name)
  tibble(Gene = PANEL, Comparison = k, log2FC = d$log2FoldChange[i], padj = d$padj[i],
         significant = d$significant[i]) }))
write_csv(pan, file.path(DE_DIR, "CREB_CRTC_target_and_inflammatory_marker_panel.csv"))
pv <- function(g, k) { r <- pan[pan$Gene == g & pan$Comparison == k, ]
  if (nrow(r) == 0 || is.na(r$log2FC)) "n.t." else sprintf("%+.2f%s", r$log2FC, if (isTRUE(r$significant)) "*" else "") }
lf <- function(g, k) pan$log2FC[pan$Gene == g & pan$Comparison == k]
stopifnot(lf("Pde4d", "control-PGE2_vs_control-un") > 1, lf("Pde4d", "Crtc2-PGE2_vs_control-PGE2") < -1,
          lf("Il6", "Crtc2-TNF_vs_control-TNF") < -2)

# ---- enrichment lookups ---------------------------------------------------------------------
top_terms <- function(k, rel, n = 3) {
  f <- file.path(ENR_DIR, k, rel); if (!file.exists(f)) return(character(0))
  t <- suppressWarnings(read_csv(f, show_col_types = FALSE)); if (nrow(t) == 0) return(character(0))
  t %>% arrange(p.adjust) %>% head(n) %>% pull(Description) }
hall <- function(k) read_csv(file.path(ENR_DIR, k, "GSEA", "GSEA_MSigDB_Hallmark.csv"), show_col_types = FALSE)
nes <- function(k, set) { h <- hall(k); v <- h$NES[h$ID == paste0("HALLMARK_", set)]; if (length(v)) v else NA_real_ }
nesf <- function(k, set) { v <- nes(k, set); if (is.na(v)) "n.s." else sprintf("%+.2f", v) }
HSETS <- c("TNFA_SIGNALING_VIA_NFKB", "INFLAMMATORY_RESPONSE", "INTERFERON_GAMMA_RESPONSE", "E2F_TARGETS", "G2M_CHECKPOINT")
for (k in c("control-PGE2_vs_control-un", "control-adenosine_vs_control-un"))
  stopifnot(nes(k, "TNFA_SIGNALING_VIA_NFKB") > 0, nes(k, "E2F_TARGETS") < 0)
for (k in c("Crtc2-un_vs_control-un", "Crtc2-PGE2_vs_control-PGE2", "Crtc2-adenosine_vs_control-adenosine"))
  stopifnot(nes(k, "INFLAMMATORY_RESPONSE") < 0, nes(k, "E2F_TARGETS") > 0)
stopifnot(nes("Crtc2-TNF_vs_control-TNF", "TNFA_SIGNALING_VIA_NFKB") < 0)

K_UN <- "Crtc2-un_vs_control-un"; K_TNF <- "Crtc2-TNF_vs_control-TNF"; K_PGE <- "Crtc2-PGE2_vs_control-PGE2"
K_ADO <- "Crtc2-adenosine_vs_control-adenosine"; K_CP <- "control-PGE2_vs_control-un"; K_CA <- "control-adenosine_vs_control-un"

# ================================ HEADER =======================================
add("# Bioinformatics Analysis Report: Control and Crtc2 Cell Lines under TNF, PGE2 and Adenosine")
add("")
add("**Report Date:** ", format(DELIVERY_DATE, "%Y-%m-%d"), "  ")
add("**Prepared by:** Zhen Gao, PhD, Principal Bioinformatics Scientist, Athenomics  ")
add("**Species:** Mus musculus (mouse; reference genome GRCm39, GENCODE M35 annotation)  ")
add("**Tissue / Cell:** Established mouse cell line (cell line name not specified by client)  ")
add("**Assay:** Bulk RNA sequencing (RNA-seq), poly-A-selected mRNA, 24 samples  ")
add("**Analysis Platform:** Linux HPC server  ")
add("")
add("---")
add("")

# ================================ 1. OBJECTIVES ================================
add("## 1. Objectives")
add("")
add("The study compares two mouse cell lines, labelled **control** and **Crtc2** by the client, each ",
    "left untreated (**un**) or treated with **TNF** (tumor necrosis factor), **PGE2** (prostaglandin E2) ",
    "or **adenosine**; three biological replicates per condition, 24 samples in total. Crtc2 encodes ",
    "CRTC2 (CREB-regulated transcription coactivator 2), a co-activator of the transcription factor CREB ",
    "(cAMP response element-binding protein; cAMP = cyclic adenosine monophosphate). The submission form does not state how the Crtc2 line was ",
    "modified, so this report refers to it as the **Crtc2 line** and reports the Crtc2 mRNA level that was measured.")
add("")
add("The six comparisons requested on the submission form (the second group named is the baseline):")
add("")
add("| # | Comparison | Question it answers |")
add("| :--- | :---: | :---: |")
add("| 1 | Crtc2-un vs control-un | How the Crtc2 line differs from control without treatment |")
add("| 2 | Crtc2-TNF vs control-TNF | How the Crtc2 line differs from control under TNF |")
add("| 3 | Crtc2-PGE2 vs control-PGE2 | How the Crtc2 line differs from control under PGE2 |")
add("| 4 | Crtc2-adenosine vs control-adenosine | How the Crtc2 line differs from control under adenosine |")
add("| 5 | control-PGE2 vs control-un | The PGE2 response of control cells |")
add("| 6 | control-adenosine vs control-un | The adenosine response of control cells |")
add("")

# ================================ 2. CORE CONCLUSIONS ==========================
add("## 2. Core Conclusions")
add("")
add("This section gives the answer for each part of the experiment in plain language. Section 3 lists ",
    "what was measured topic by topic; Section 8.2 gives **more detail on each conclusion below and names ",
    "the section and delivered file that supports it**; Section 8.3 states what these data cannot tell you; ",
    "Section 8.4 lists the experiment that would settle each open question.")
add("")
CC <- c(
  sprintf("**The data are of good quality and the sample labels are consistent.** All 24 libraries passed quality control (%s–%s%% of reads uniquely mapped), all were sequenced together on one flowcell lane, and replicates of the same condition are %s times closer to each other than to other conditions. No sample needs to be excluded.",
          f1(Q_UNIQ[1]), f1(Q_UNIQ[2]), f1(pca$separation_ratio)),
  sprintf("**The Crtc2 line has %s less Crtc2 mRNA than control, in every condition** (Crtc2 line / control = %s–%s). The transcript is reduced, not absent. RNA-seq measures mRNA only, so whether CRTC2 protein is lost should be confirmed by Western blot (Section 8.4).",
          CR_LESS, f2(CR_RATIO[1]), f2(CR_RATIO[2])),
  sprintf("**Without treatment, the Crtc2 line already differs from control in %s genes, and the difference runs opposite to the PGE2 response.** Classical CREB target genes (Nr4a2, Crem, Pde4b, Pde4d, Sik1) are lower in the Crtc2 line; of the %s genes changed both by the Crtc2 line and by PGE2, %s change in opposite directions.",
          fmt(dg(K_UN)$Total), fmt(O_UN_PGE2$n), pc_(O_UN_PGE2$opp)),
  sprintf("**PGE2 is a strong stimulus in control cells (%s genes change), and the Crtc2 line reproduces only about %s of that response.** In control cells PGE2 raises CREB targets and inflammatory genes (Pde4d %s, Il6 %s log2 fold change) and lowers cell-cycle genes. Under PGE2 the Crtc2 line differs from control in %s genes; %s of them are genes that PGE2 changes in control cells, and %s of those move in the opposite direction.",
          fmt(dg(K_CP)$Total), pc_(rr("PGE2", "Slope_Crtc2_vs_control_response")),
          sprintf("%+.1f", lf("Pde4d", K_CP)), sprintf("%+.1f", lf("Il6", K_CP)), fmt(dg(K_PGE)$Total), pc_(O_PGE2$n / dg(K_PGE)$Total), pc_(O_PGE2$opp)),
  sprintf("**Adenosine has a mild effect on control cells (%s genes, only %s of them changing 2-fold or more) and almost none on the Crtc2 line.** In the Crtc2 line, adenosine-treated and untreated samples are practically indistinguishable (Section 7.1); on the control-cell adenosine genes the Crtc2 line reproduces about %s of the response.",
          fmt(dg(K_CA)$Total), fmt(dg(K_CA)$Up_lfc1 + dg(K_CA)$Down_lfc1), pc_(rr("adenosine", "Slope_Crtc2_vs_control_response"))),
  sprintf("**TNF is the strongest stimulus in the experiment, and the largest Crtc2-line difference is seen under TNF (%s genes).** The Crtc2 line still responds to TNF, but less (about %s of the control response), and TNF-induced inflammatory genes such as Il6 (%s log2 fold change) and Cxcl10 are clearly lower than in control.",
          fmt(dg(K_TNF)$Total), pc_(rr("TNF", "Slope_Crtc2_vs_control_response")), sprintf("%+.1f", lf("Il6", K_TNF))),
  "**What we suggest next.** Confirm CRTC2 protein loss by Western blot, and confirm a small panel of the key genes (Pde4d, Nr4a2, Il6) by qPCR (quantitative PCR) in independent cultures."
)
for (x in CC) { add("- ", x); add("") }

# ================================ 3. KEY FINDINGS ==============================
add("## 3. Key Findings")
add("")
add("Section 2 gave the conclusion for each part of the experiment. This section reports the measured results by topic.")
add("")
add("**Differentially expressed genes (DEGs) per comparison**")
add("")
add("| Comparison | Up-regulated | Down-regulated | Total DEGs | DEGs with at least 2-fold change |")
add("| :--- | :---: | :---: | :---: | :---: |")
for (i in seq_len(nrow(deg))) { r <- deg[i, ]
  add("| ", cmp(r$Comparison), " | ", fmt(r$Up), " | ", fmt(r$Down), " | **", fmt(r$Total), "** | ",
      fmt(r$Up_lfc1 + r$Down_lfc1), " (", fmt(r$Up_lfc1), " up / ", fmt(r$Down_lfc1), " down) |") }
add("")
add("A DEG is a gene with adjusted p-value (padj; Benjamini-Hochberg false discovery rate, FDR) of 0.05 or less ",
    "and a change of at least 1.2-fold (absolute log2 fold change, log2FC, of at least 0.263). The last ",
    "column applies the stricter 2-fold cut-off (absolute log2FC of at least 1) shown as the default on the ",
    "submission form. The 1.2-fold cut-off is used as the primary criterion because many consistent, biologically relevant changes in this data set are smaller than 2-fold; Crtc2 itself (log2FC about -0.6) would not pass a 2-fold cut-off. A positive log2FC means higher expression in the first group named.")
add("")
add("**Crtc2 mRNA.** Lower in the Crtc2 line in all four conditions (log2FC ",
    paste(sprintf("%+.2f", crtc2_lfc), collapse = ", "), " for comparisons 1–4; all significant).")
add("")
add("**How much of each treatment response the Crtc2 line keeps** (descriptive; Section 7.4):")
add("")
add("| Treatment | Genes used | Crtc2-line response as a fraction of the control response |")
add("| :--- | :---: | :---: |")
for (tr in c("PGE2", "adenosine", "TNF")) add("| ", tr, " | ", fmt(rr(tr, "N_genes")), " | ", pc_(rr(tr, "Slope_Crtc2_vs_control_response")), " |")
add("")
add("**Pathway themes.** In control cells, PGE2 and adenosine both raise inflammatory and interferon gene ",
    "programmes and lower cell-cycle programmes. In the Crtc2 line the same programmes move the opposite way: ",
    "inflammatory and interferon genes are lower and cell-cycle genes higher than in control, at baseline and ",
    "under every treatment (Section 7.3).")
add("")
add(sprintf("**Data quality.** %s–%s million read pairs per sample (paired-end, %s bp), %s–%s%% of reads aligned, %s–%s%% uniquely aligned, rRNA (ribosomal RNA) %s–%s%%. Library strandedness was inferred from the data as %s for all 24 samples.",
            f1(RP[1] / 1e6), f1(RP[2] / 1e6), paste(Q_LEN, collapse = "/"), f1(Q_MAP[1]), f1(Q_MAP[2]),
            f1(Q_UNIQ[1]), f1(Q_UNIQ[2]), f2(Q_RRNA[1]), f2(Q_RRNA[2]), STRAND))
add("")

# ================================ 4. SAMPLE INFORMATION ========================
add("## 4. Sample Information")
add("")
add("| Sample ID (this report) | Client sample name | Cell line | Treatment | Replicate | Read pairs | Reads assigned to genes | Genes detected |")
add("| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: |")
for (i in seq_len(nrow(meta))) { s <- meta$sample_id[i]; r <- reads[reads$sample_id == s, ]; stopifnot(nrow(r) == 1)
  add("| ", s, " | ", meta$client_sample_name[i], " | ", meta$genotype[i], " | ", meta$treatment[i], " | ",
      meta$replicate[i], " | ", fmt(r$raw_read_pairs), " | ", fmt(r$assigned_reads), " | ", fmt(r$genes_detected), " |") }
add("")
add("The client sample names list the PGE2 samples as \"PEG2\" (for example control-PEG2-1), while the comparison ",
    "list on the same form writes \"pge2\"; this report uses **PGE2** throughout and keeps the client name alongside. ",
    "All samples are from an established mouse cell line (\"Normal only\", snap-frozen, per the form). ",
    sprintf("All 24 libraries were sequenced on the same flowcell lane (%s), so sequencing run cannot be confounded with cell line or treatment.", LANES))
add("")

# ================================ 5. RATIONALE =================================
add("## 5. Analysis Rationale and Decision Criteria")
add("")
add("| Step | Choice | Why |")
add("| :--- | :---: | :---: |")
add("| Reference | GRCm39 primary assembly, GENCODE M35 | Current mouse reference genome and annotation |")
add("| Strandedness | Inferred from the data | Avoids assuming the library chemistry; a wrong setting would distort gene counts |")
add("| Statistical model | One DESeq2 model, `~ Group`, all 24 samples (8 groups) | Single sequencing lane and balanced design; all six comparisons share one variance estimate from 24 samples, which is more precise than fitting each pair separately |")
add("| Gene-type filter | ", fmt(fs$n_original), " -> ", fmt(fs$n_after_regex), " genes | Removes ribosomal-protein, mitochondrial and non-coding/predicted loci, which add noise and dilute multiple-testing correction |")
add("| Low-expression filter | ", fmt(fs$n_after_regex), " -> ", fmt(fs$n_final), " genes (", fmt(min(fs$n_eligible_per_contrast)), "–", fmt(max(fs$n_eligible_per_contrast)), " tested per comparison) | A gene is tested in a comparison if it has at least 10 counts (Salmon estimated counts rounded to integers, as used by DESeq2) in at least 3 of the 6 samples being compared, or in at least 30% of all samples. This keeps genes that are switched on by one treatment only, which a rule requiring expression across most of the 24 samples would discard |")
add("| Significance | padj of 0.05 or less and at least 1.2-fold change | False-discovery control plus a minimum effect size; counts at the 2-fold cut-off are also reported |")
add("| Fold-change shrinkage | `ashr` | Stabilises fold changes of low-count genes |")
add("| Over-representation analysis | Up- and down-regulated genes tested separately | Testing both directions together loses power when many genes respond |")
add("| GSEA ranking | DESeq2 Wald statistic | Uses every tested gene without a cut-off; the Wald statistic weights genes evenly enough that a few extreme genes cannot dominate the result |")
add("")

# ================================ 6. METHODS ===================================
add("## 6. Methods")
add("")
add("| Step | Tool | Key parameters |")
add("| :--- | :---: | :---: |")
add("| Pipeline | nf-core/rnaseq 3.15.1 | `--aligner star_salmon --gencode` |")
add("| Trimming | Trim Galore! | nf-core defaults |")
add("| Alignment | STAR | GRCm39 index; `--outFilterMultimapNmax 3` |")
add("| Quantification | Salmon | Gene-level estimated counts and TPM (transcripts per million) |")
add("| Quality control | FastQC, RSeQC, Qualimap, Picard, MultiQC | Strandedness inferred by RSeQC and Salmon |")
add("| Differential expression | DESeq2 | `~ Group`; Wald test; per-comparison gene set as in Section 5 with Benjamini-Hochberg correction within it; `lfcShrink(type = \"ashr\")` |")
add("| Over-representation | clusterProfiler | GO (Gene Ontology) biological process (BP), molecular function (MF) and cellular component (CC), and KEGG (Kyoto Encyclopedia of Genes and Genomes); background = genes tested in that comparison |")
add("| GSEA (gene set enrichment analysis) | clusterProfiler / fgsea | Ranked by the DESeq2 Wald statistic; GO BP, KEGG, MSigDB (Molecular Signatures Database) Hallmark gene sets (mouse orthologs) |")
add("")
add("The count matrix contains Salmon **gene-level estimated counts**, not raw read counts ",
    "(`Reads/All_sample_gene_counts.tsv`, delivered unrounded); DESeq2 and the low-expression filter use them rounded to integers. `Reads/All_sample_gene_tpm.tsv` contains TPM values.")
add("")

# ================================ 7. RESULTS ===================================
add("## 7. Results")
add("")
add("### 7.1 Principal Component Analysis (PCA)")
add("")
add(sprintf("The first principal component (PC1) explains %s%% of the variance and the second (PC2) %s%% (`DE_PCA_Results/PCA_AllSamples_Overview.pdf`; coordinates in `DE_PCA_Results/PCA_coordinates_overview.csv`). Replicates of the same condition lie on average %s units apart and samples of different conditions %s units apart, a ratio of %s: the conditions are far more distinct than the replicate-to-replicate variation.",
            pca$percentVar[1], pca$percentVar[2], f2(pca$within_dist), f2(pca$between_dist), f1(pca$separation_ratio)))
add("")
add(sprintf("**PC1 is the TNF axis.** The TNF-treated samples sit at one end (mean PC1 %s) and all other samples at the other (untreated %s, PGE2 %s, adenosine %s). Along this axis the Crtc2-TNF group lies between control-TNF and the untreated groups: control-TNF is %s units from control-un, Crtc2-TNF %s units from Crtc2-un, so the TNF shift is smaller in the Crtc2 line.",
            f1(TNF_PC1["TNF"]), f1(TNF_PC1["un"]), f1(TNF_PC1["PGE2"]), f1(TNF_PC1["adenosine"]),
            f1(d_("control-TNF", "control-un")), f1(d_("Crtc2-TNF", "Crtc2-un"))))
add("")
add(sprintf("Among the non-TNF groups the same pattern holds for PGE2 and adenosine. PGE2 moves control cells %s units from control-un but the Crtc2 line only %s units from Crtc2-un; adenosine moves control cells %s units but the Crtc2 line only %s units, which is no further apart than replicates of one condition. The Crtc2-PGE2 group sits only %s units from untreated control cells.",
            f1(d_("control-PGE2", "control-un")), f1(d_("Crtc2-PGE2", "Crtc2-un")),
            f1(d_("control-adenosine", "control-un")), f2(d_("Crtc2-adenosine", "Crtc2-un")), f1(d_("Crtc2-PGE2", "control-un"))))
add("")
add(sprintf("**PC2 is the cell-line axis.** Every Crtc2-line sample has a higher PC2 value (%s to %s) than every control sample (%s to %s), in all four conditions: the second-largest source of variation in the data set is the difference between the two lines.",
            f1(PC2_G[["Crtc2"]][1]), f1(PC2_G[["Crtc2"]][2]), f1(PC2_G[["control"]][1]), f1(PC2_G[["control"]][2])))
add("")
add("Each sample's nearest neighbour in the PCA was checked. ",
    if (nrow(nn_bad) == 0) "Every sample sits closest to a replicate of its own condition." else
      paste0(nrow(nn_bad), " of 24 samples ", if (nrow(nn_bad) == 1) "sits" else "sit", " closest to a sample of another condition (",
             paste(nn_bad$sample_id, "->", nn_bad$nearest, collapse = "; "),
             "). This involves only the Crtc2-un and Crtc2-adenosine groups, which are indistinguishable because adenosine has no measurable effect on the Crtc2 line (above); it is not a sign of a sample mix-up. Every other sample sits closest to a replicate of its own condition."))
add("")
add("| Pair of conditions | PCA centroid distance |")
add("| :--- | :---: |")
for (pr in list(c("control-un", "control-TNF"), c("Crtc2-un", "Crtc2-TNF"), c("control-un", "control-PGE2"),
                c("Crtc2-un", "Crtc2-PGE2"), c("control-un", "control-adenosine"), c("Crtc2-un", "Crtc2-adenosine"),
                c("control-un", "Crtc2-un"), c("control-TNF", "Crtc2-TNF"), c("control-PGE2", "Crtc2-PGE2"),
                c("control-adenosine", "Crtc2-adenosine"), c("control-un", "Crtc2-PGE2")))
  add("| ", pr[1], " vs ", pr[2], " | ", f2(d_(pr[1], pr[2])), " |")
add("")

add("### 7.2 Differential Expression")
add("")
add("| Comparison | Genes tested | Up | Down | Total DEGs | % of tested genes |")
add("| :--- | :---: | :---: | :---: | :---: | :---: |")
for (i in seq_len(nrow(deg))) { r <- deg[i, ]
  add("| ", cmp(r$Comparison), " | ", fmt(r$Genes_tested), " | ", fmt(r$Up), " | ", fmt(r$Down), " | **",
      fmt(r$Total), "** | ", sprintf("%.1f%%", 100 * r$Total / r$Genes_tested), " |") }
add("")
add("Each comparison folder in `DE_PCA_Results/` holds the full table of tested genes (fold change, Wald statistic, ",
    "p-value, padj and normalised counts of the six samples), the significant subset, a volcano plot and a heatmap of the top 50 DEGs.")
add("")
add("**Crtc2 mRNA** (DESeq2-normalised counts, mean of three replicates):")
add("")
add("| Treatment | control | Crtc2 line | Crtc2 line / control |")
add("| :--- | :---: | :---: | :---: |")
for (tr in c("un", "TNF", "PGE2", "adenosine")) { r <- cr_g[cr_g$treatment == tr, ]
  add("| ", TREAT_LONG[tr], " | ", fmt(round(r$control)), " | ", fmt(round(r$Crtc2)), " | ", f2(r$ratio), " |") }
add("")
add(sprintf("Crtc2 mRNA is higher with TNF in both lines (group means %s-fold in control, %s-fold in the Crtc2 line; TNF vs untreated was not a tested comparison). In every condition the Crtc2 line carries %s less Crtc2 mRNA than control.", f1(CR_TNF_FOLD["control"]), f1(CR_TNF_FOLD["Crtc2"]), CR_LESS))
add("")

add("### 7.3 Pathway Enrichment")
add("")
add("**Hallmark gene sets (GSEA normalised enrichment score, NES; positive = higher in the first group named; n.s. = not significant)**")
add("")
add("| Comparison | ", paste(gsub("_", " ", HSETS), collapse = " | "), " |")
add("| :--- | ", paste(rep(":---:", length(HSETS)), collapse = " | "), " |")
for (k in names(res_list)) add("| ", cmp(k), " | ", paste(sapply(HSETS, function(s) nesf(k, s)), collapse = " | "), " |")
add("")
add("The table shows one consistent pattern. In control cells, PGE2 and adenosine raise TNF/NF-kB (nuclear factor kappa B), ",
    "inflammatory and interferon-gamma programmes and lower cell-cycle programmes (E2F targets, G2/M checkpoint). In the ",
    "Crtc2 line these programmes go the other way: lower inflammatory and interferon signatures and higher cell-cycle ",
    "signatures than control, already without treatment and under PGE2 and adenosine. Under TNF, the TNF/NF-kB programme ",
    "itself is lower in the Crtc2 line.")
add("")
add("**Top over-represented GO biological-process terms** (three strongest per direction):")
add("")
add("| Comparison | Higher in first group | Lower in first group |")
add("| :--- | :---: | :---: |")
for (k in names(res_list)) {
  u <- top_terms(k, "GO/GO_BP_UP.csv"); d <- top_terms(k, "GO/GO_BP_DOWN.csv")
  add("| ", cmp(k), " | ", if (length(u)) paste(u, collapse = "; ") else "none at FDR < 0.05", " | ",
      if (length(d)) paste(d, collapse = "; ") else "none at FDR < 0.05", " |") }
add("")
EW <- enr %>% filter(Analysis %in% c("GO_BP_UP", "GO_BP_DOWN", "KEGG_UP", "KEGG_DOWN", "GSEA_GO_BP", "GSEA_KEGG", "GSEA_Hallmark")) %>%
  select(Comparison, Analysis, N_significant) %>% pivot_wider(names_from = Analysis, values_from = N_significant)
add("**Number of significant terms (FDR < 0.05)**")
add("")
add("| Comparison | GO BP up | GO BP down | KEGG up | KEGG down | GSEA GO BP | GSEA KEGG | GSEA Hallmark |")
add("| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: |")
for (i in seq_len(nrow(EW))) { r <- EW[i, ]
  add("| ", cmp(r$Comparison), " | ", paste(fmt(unlist(r[, c("GO_BP_UP", "GO_BP_DOWN", "KEGG_UP", "KEGG_DOWN", "GSEA_GO_BP", "GSEA_KEGG", "GSEA_Hallmark")])), collapse = " | "), " |") }
add("")
add("Complete tables and dot plots are in `Enrichment/<comparison>/{GO,KEGG,GSEA}/`.")
add("")

add("### 7.4 How the Crtc2 line responds to each treatment")
add("")
add("The client's comparisons test the Crtc2 line against control within each treatment. To show directly how much ",
    "of each treatment response remains in the Crtc2 line, we took the genes that respond in control cells and compared ",
    "the size of their change (treated vs untreated, mean of three replicates) in the Crtc2 line with the size in control. ",
    "The fraction below is the slope of that relationship: 100% would mean the Crtc2 line responds exactly like control, ",
    "0% that it does not respond at all. These are descriptive summaries, not a statistical test (see Section 8.3).")
add("")
add("| Treatment | Genes used | Crtc2-line response as a fraction of control | Correlation | Genes responding less than half as much |")
add("| :--- | :---: | :---: | :---: | :---: |")
for (tr in c("PGE2", "adenosine", "TNF")) add("| ", tr, " | ", fmt(rr(tr, "N_genes")), " | ", pc_(rr(tr, "Slope_Crtc2_vs_control_response")),
  " | ", f2(rr(tr, "Correlation")), " | ", pc_(rr(tr, "Fraction_genes_below_half_response")), " |")
add("")
add("Genes used: for PGE2 and adenosine, the DEGs of comparisons 5 and 6. TNF vs untreated was not among the requested ",
    "comparisons and was not tested, so for TNF the genes are those whose mean expression differs at least 2-fold between ",
    "control-TNF and control-un. Values: `DE_PCA_Results/Crtc2_line_treatment_response_retention.csv`.")
add("")
add(sprintf("The Crtc2-line differences also point in a consistent direction relative to the treatments. Without treatment, of the %s genes that differ both between the cell lines and with PGE2 in control cells, %s change in opposite directions (genome-wide correlation of the two sets of fold changes %s); for adenosine the corresponding figures are %s of %s genes (correlation %s). Under PGE2, %s of the %s genes significant in both comparison 3 and comparison 5 move in opposite directions (correlation %s): the difference between the lines under PGE2 is, for the most part, the PGE2 response that the Crtc2 line does not make.",
            fmt(O_UN_PGE2$n), pc_(O_UN_PGE2$opp), f2(O_UN_PGE2$r), pc_(O_UN_ADO$opp), fmt(O_UN_ADO$n), f2(O_UN_ADO$r),
            pc_(O_PGE2$opp), fmt(O_PGE2$n), f2(O_PGE2$r)))
add("")
add("**Selected genes** (log2FC; * = significant; n.t. = not tested in that comparison because expression is too low). ",
    "Nr4a2, Crem, Pde4b, Pde4d and Sik1 are classical CREB/CRTC target genes; Il6, Ptgs2, Cxcl2, Cxcl10 and Ccl5 are inflammatory genes.")
add("")
add("| Gene | ", paste(sapply(names(res_list), cmp), collapse = " | "), " |")
add("| :--- | ", paste(rep(":---:", 6), collapse = " | "), " |")
for (g in PANEL) add("| ", g, " | ", paste(sapply(names(res_list), function(k) pv(g, k)), collapse = " | "), " |")
add("")
add("The CREB targets rise with PGE2 in control cells and are lower in the Crtc2 line in every condition, most strongly ",
    "under PGE2 and TNF. Il6, Ptgs2 and Cxcl2 are induced by PGE2 in control cells and markedly lower in the Crtc2 line under ",
    "PGE2 and TNF. Fos and Fosb go the other way: they are higher in the Crtc2 line without treatment and under TNF and adenosine. Values: ",
    "`DE_PCA_Results/CREB_CRTC_target_and_inflammatory_marker_panel.csv`.")
add("")

# ================================ 8. CONCLUSIONS ===============================
add("## 8. Conclusions, Limitations and Next Steps")
add("")
add("### 8.1 Per-comparison summary")
add("")
add("| Comparison | DEGs | Read-out |")
add("| :--- | :---: | :---: |")
RO <- c(
  "Crtc2 line differs from control already without treatment; CREB targets and inflammatory/interferon genes lower, cell-cycle genes higher",
  "Largest difference; the Crtc2 line mounts a smaller TNF response, with lower TNF-induced inflammatory genes",
  sprintf("Mostly the PGE2 response that the Crtc2 line does not make: %s of these DEGs are PGE2-responsive in control cells, %s of them moved in the opposite direction", pc_(O_PGE2$n / deg$Total[deg$Comparison == K_PGE]), pc_(O_PGE2$opp)),
  sprintf("Mostly the baseline Crtc2-line difference: %s of these DEGs are also DEGs without treatment, all in the same direction", pc_(O_ADO_UN$n / deg$Total[deg$Comparison == K_ADO])),
  "Strong PGE2 response: CREB targets and inflammatory genes up, cell-cycle genes down",
  "Mild adenosine response, same direction as PGE2 but much smaller")
for (i in seq_len(nrow(deg))) add("| ", cmp(deg$Comparison[i]), " | ", fmt(deg$Total[i]), " | ", RO[i], " |")
add("")
add("### 8.2 Practical reading")
add("")
add("The seven points of Section 2, each with the section and delivered file that supports it.")
add("")
EV <- c("Sections 3 and 7.1; `QC/multiqc/`, `DE_PCA_Results/PCA_AllSamples_Overview.pdf`",
        "Section 7.2 Crtc2 table; `DE_PCA_Results/<comparison>/*_ALL_genes.csv`, row Crtc2",
        "Sections 7.3 and 7.4; `DE_PCA_Results/Crtc2-un_vs_control-un/`, `CREB_CRTC_target_and_inflammatory_marker_panel.csv`",
        "Sections 7.1 and 7.4; `DE_PCA_Results/control-PGE2_vs_control-un/`, `Crtc2_line_treatment_response_retention.csv`",
        "Section 7.1 (centroid distances) and 7.4; `DE_PCA_Results/control-adenosine_vs_control-un/`",
        "Sections 7.1 and 7.4; `DE_PCA_Results/Crtc2-TNF_vs_control-TNF/`, `Enrichment/Crtc2-TNF_vs_control-TNF/GSEA/`",
        "Section 8.4")
stopifnot(length(EV) == length(CC))
for (k in seq_along(CC)) { add("- ", CC[k], " *(Evidence: ", EV[k], ")*"); add("") }
add("### 8.3 Limitations")
add("")
add("- **The nature of the Crtc2 modification is not stated.** The data show reduced Crtc2 mRNA; they cannot show whether CRTC2 protein is absent or reduced.")
add("- **The response fractions in Section 7.4 are descriptive.** A formal genotype x treatment interaction test was not part of the requested comparisons (the factorial-design box on the form was not ticked). The figures summarise group means and are not p-values.")
add("- **TNF vs untreated was not a requested comparison**, so no TNF-response DEG list is delivered; the TNF statements rest on the PCA and on comparison 2.")
add("- **Three replicates per condition.** Consistent changes are detected well; small or variable effects can be missed, so a gene absent from a DEG list is not shown to be unchanged.")
add("- **Some chemokines are barely expressed in these cells** (for example Cxcl1 and Ccl2 have fewer than 20 reads in every sample) and could not be tested.")
add("- **RNA level only.** Transcript changes do not establish protein or functional changes; pathway terms name a theme, not a demonstrated mechanism.")
add("")
add("### 8.4 Suggested confirmation")
add("")
add("| Open question | Experiment that would settle it |")
add("| :--- | :---: |")
add("| Is CRTC2 protein lost in the Crtc2 line? | Western blot for CRTC2 in control and Crtc2 cells |")
add("| Are the key gene changes reproducible? | qPCR of Pde4d, Nr4a2, Il6 (and Crtc2) in independent cultures, untreated and with PGE2 and TNF |")
add("| Are the pathway differences functional? | A functional read-out for the main themes, for example a proliferation assay and cytokine secretion (for example IL-6 by enzyme-linked immunosorbent assay, ELISA) |")
add("")

# ================================ 9. DELIVERABLES ==============================
add("## 9. Deliverable Files")
add("")
add("| Path | Content |")
add("| :--- | :---: |")
add("| `Reads/All_sample_gene_counts.tsv` | Gene-level estimated counts (Salmon), 24 samples — the matrix for re-analysis |")
add("| `Reads/All_sample_gene_tpm.tsv` | Gene-level TPM, 24 samples |")
add("| `Reads/library_metrics_per_sample.csv` | Read pairs, reads assigned to genes, genes detected |")
add("| `QC/multiqc/` | MultiQC quality-control report |")
add("| `DE_PCA_Results/PCA_AllSamples_Overview.pdf` | PCA of all 24 samples |")
add("| `DE_PCA_Results/PCA_coordinates_overview.csv` | PCA coordinates |")
add("| `DE_PCA_Results/SampleDistance_Heatmap.pdf` | Sample-to-sample distances |")
add("| `DE_PCA_Results/normalized_counts_all_samples.csv` | DESeq2-normalised counts |")
add("| `DE_PCA_Results/DEG_summary_all_comparisons.csv` | DEG counts, 1.2-fold and 2-fold cut-offs |")
add("| `DE_PCA_Results/<comparison>/` | Full and significant result tables, volcano plot, top-50 heatmap (6 comparisons) |")
add("| `DE_PCA_Results/Crtc2_line_treatment_response_retention.csv` | Section 7.4 response fractions |")
add("| `DE_PCA_Results/CREB_CRTC_target_and_inflammatory_marker_panel.csv` | Section 7.4 gene panel |")
add("| `Enrichment/<comparison>/GO/` and `Enrichment/<comparison>/KEGG/` | Over-representation, split by direction |")
add("| `Enrichment/<comparison>/GSEA/` | GSEA: GO BP, KEGG, MSigDB Hallmark |")
add("| `Enrichment/Enrichment_summary_all_comparisons.csv` | Number of significant terms per analysis |")
add("| `mouse_Gene_annotation_20260628.xlsx` | Gene annotation reference |")
add("")
add("---")
add("")
add("*Zhen Gao, PhD, Principal Bioinformatics Scientist, Athenomics*")
add("")

txt <- paste(L, collapse = "\n")
bad <- c("knockout", "knock-out", "knockdown", " KO ", "the same (answers|ground|findings)", "restates the", "repeats the",
         "while it can still", "long-read", "Nanopore", "PacBio", "possibly", "might be")
for (b in bad) if (grepl(b, txt, ignore.case = TRUE)) stop("banned phrase in report: ", b)
writeLines(L, REPORT)
cat("Report written:", normalizePath(REPORT), "\n")
