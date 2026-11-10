#!/usr/bin/env Rscript
# P52 step 6a: read existing results and export report facts; no statistical refit.
suppressPackageStartupMessages({library(jsonlite);library(ggplot2);library(ggrepel);library(DESeq2)})
args <- commandArgs(trailingOnly=FALSE)
S <- dirname(normalizePath(sub('^--file=','',grep('^--file=',args,value=TRUE)[1])))
setwd(S); D <- '../Data_Analysis_20261109'; C <- 'rds_cache'
out <- '../logs/report_facts_20261109.json'; stopifnot(!file.exists(out))
pca <- readRDS(file.path(C,'pca_summary.rds')); meta <- readRDS(file.path(C,'meta.rds')); fs <- readRDS(file.path(C,'filter_stats.rds'))
cols <- readRDS(file.path(C,'group_colors.rds'))
stopifnot(nrow(pca$pca_data)==24, setequal(pca$pca_data$name,meta$sample_id), all(table(meta$Group)==3))
p <- ggplot(pca$pca_data,aes(PC1,PC2,color=Group,shape=genotype))+geom_point(size=4)+geom_text_repel(aes(label=name),max.overlaps=Inf,seed=42)+scale_color_manual(values=cols)+scale_shape_manual(values=c(control=16,Crtc2=17))
b <- ggplot_build(p)
stopifnot(nrow(b$data[[1]])==24,all(is.finite(b$data[[1]]$x)),all(is.finite(b$data[[1]]$y)),!anyNA(b$data[[1]]$shape),!anyNA(b$data[[1]]$colour),nrow(b$data[[2]])==24,setequal(b$data[[2]]$label,meta$sample_id))
dds <- readRDS(file.path(C,'dds.rds')); nc <- counts(dds,normalized=TRUE)
stopifnot(setequal(colnames(nc),meta$sample_id))
# The source figure uses match() for both annotations; assert all labels against SSOT.
ann <- as.character(meta$Group[match(colnames(nc),meta$sample_id)])
stopifnot(!anyNA(ann),identical(ann,unname(setNames(as.character(meta$Group),meta$sample_id)[colnames(nc)])))
qc <- read.delim(file.path(D,'QC/multiqc/star_salmon/multiqc_report_data/multiqc_general_stats.txt'),check.names=FALSE)
getrange <- function(pat){x<-grep(pat,names(qc),value=TRUE);stopifnot(length(x)==1);v<-as.numeric(qc[[x]]);v<-v[is.finite(v)];stopifnot(length(v)==24);range(v)}
strand <- read.delim(file.path(D,'QC/multiqc/star_salmon/multiqc_report_data/multiqc_fail_strand_check_table.txt'),check.names=FALSE)
facts <- list(pca=list(percentVar=pca$percentVar,within_dist=pca$within_dist,between_dist=pca$between_dist,separation_ratio=pca$separation_ratio,nearest=pca$nn_tbl,plotted_points=nrow(b$data[[1]]),plotted_labels=nrow(b$data[[2]])),filters=fs,qc=list(mapped=getrange('star_mapped_percent_1$'),unique=getrange('star_uniquely_mapped_percent_1$'),rrna=getrange('percent_rRNA$'),strand=unique(strand[['Inferred strandedness']])),versions=sapply(c('DESeq2','ashr','clusterProfiler','fgsea','org.Mm.eg.db'),function(x)as.character(packageVersion(x))),annotation_check='24 labels match sample SSOT; PCA 24 finite points and 24 text labels; heatmap match-based annotation verified')
write_json(facts,out,auto_unbox=TRUE,pretty=TRUE,digits=15)
cat('P52 report facts exported; 24/24 points and labels; sample annotation assertions PASS\n')
