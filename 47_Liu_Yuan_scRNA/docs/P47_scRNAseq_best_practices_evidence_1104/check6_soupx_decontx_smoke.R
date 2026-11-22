suppressPackageStartupMessages({library(SoupX); library(celda); library(SingleCellExperiment); library(Matrix)})
set.seed(1); g <- 300; n <- 400
prof <- matrix(rgamma(g*3, 0.5), g, 3); cl <- sample(1:3, n, TRUE)
soup <- rowMeans(prof); soup <- soup/sum(soup)
X <- sapply(1:n, function(i) rpois(g, 2000*(0.9*prof[,cl[i]]/sum(prof[,cl[i]]) + 0.1*soup)))
X <- Matrix(X, sparse=TRUE); rownames(X) <- paste0("g",1:g); colnames(X) <- paste0("c",1:n)
tod <- cbind(X, Matrix(rpois(g*600, soup*50), g, 600, sparse=TRUE)); colnames(tod) <- c(colnames(X), paste0("e",1:600))
sc <- SoupChannel(tod, X); sc <- setClusters(sc, cl); sc <- setContaminationFraction(sc, 0.1)
out <- adjustCounts(sc, roundToInt=TRUE)
stopifnot(all(dim(out)==dim(X)), sum(out) < sum(X)); cat("SoupX OK: removed frac =", round(1-sum(out)/sum(X),3), "\n")
sce <- SingleCellExperiment(list(counts=X)); sce <- decontX(sce, z=cl, verbose=FALSE)
cont <- colData(sce)$decontX_contamination
stopifnot(length(cont)==n, all(is.finite(cont))); cat("decontX OK: median contamination =", round(median(cont),3), "\n")
