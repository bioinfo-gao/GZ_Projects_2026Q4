suppressPackageStartupMessages({library(ClusterDE); library(countsplit); library(Matrix); library(SingleCellExperiment)})
set.seed(7); G <- 300; N <- 600
mu <- exp(rnorm(G, 0, 1)); sim <- function(lfc) {
  grp <- rep(1:2, each = N/2); M <- sapply(grp, function(k) rnbinom(G, size = 5, mu = mu * ifelse(k == 2, exp(lfc), 1)))
  rownames(M) <- paste0("g", 1:G); colnames(M) <- paste0("c", 1:N); M }
clus <- function(M) { L <- log1p(t(t(M) / colSums(M)) * 1e4); pc <- prcomp(t(L), rank. = 10)$x; kmeans(pc, 2, nstart = 10)$cluster }
pvals <- function(M, cl) { L <- log1p(t(t(M) / colSums(M)) * 1e4)
  setNames(apply(L, 1, function(x) tryCatch(t.test(x[cl == 1], x[cl == 2])$p.value, error = function(e) 1)), rownames(M)) }
run <- function(M, tag, truth = character()) {
  cl <- clus(M); p_naive <- pvals(M, cl)
  naive <- sum(p.adjust(p_naive, "BH") < 0.05)
  # count splitting: 在 train 上聚类、在 test 上检验
  cs <- countsplit(as(t(M), "CsparseMatrix"), folds = 2, overdisps = rep(5, nrow(M)))   # 要求 cell x gene
  tr <- t(as.matrix(cs[[1]])); te <- t(as.matrix(cs[[2]])); stopifnot(all(tr + te == M)); dimnames(tr) <- dimnames(te) <- dimnames(M)
  csd <- names(which(p.adjust(pvals(te, clus(tr)), "BH") < 0.05))
  # ClusterDE: 合成 null（无簇结构）上重复同一聚类+检验，做对比
  so <- suppressWarnings(Seurat::CreateSeuratObject(as(M, "CsparseMatrix")))
  nul <- constructNull(so, family = "nb", n_cores = 4, seed = 1)
  nm <- as.matrix(if (is(nul, "Seurat")) Seurat::GetAssayData(nul, layer = "counts") else if (is(nul, "SingleCellExperiment")) counts(nul) else nul); dimnames(nm) <- dimnames(M)
  res <- callDE(p_naive, pvals(nm, clus(nm)), nlogTrans = TRUE, FDR = 0.05)
  stopifnot(all(c("gene","record") %in% names(res)), nrow(res) == nrow(M)); cde <- res$gene[res$record == 1]
  hit <- function(x) if (length(truth)) sprintf(" (真 DE 命中 %d/%d)", length(intersect(x, truth)), length(truth)) else ""
  cat(sprintf("%-10s naive=%3d  countsplit=%3d%s  ClusterDE=%3d%s\n", tag, naive, length(csd), hit(csd), length(cde), hit(cde)))
}
run(sim(0), "同质(阴性)")
M2 <- sim(0); truth <- paste0("g", 1:20); M2[1:20, 301:600] <- rnbinom(20 * 300, size = 5, mu = mu[1:20] * 4)
run(M2, "20真DE(阳性)", truth)
