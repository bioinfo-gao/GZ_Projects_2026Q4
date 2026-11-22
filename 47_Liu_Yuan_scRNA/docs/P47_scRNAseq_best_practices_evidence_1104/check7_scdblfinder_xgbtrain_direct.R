suppressPackageStartupMessages(library(scDblFinder))
cat("scDblFinder", as.character(packageVersion("scDblFinder")), " xgboost", as.character(packageVersion("xgboost")), "\n")
set.seed(1); n <- 2000; X <- matrix(rnorm(n * 10), n); y <- factor(ifelse(X[,1] + rnorm(n, 0, .5) > 0, "doublet", "real"), levels = c("real","doublet"))
r <- tryCatch({ fit <- scDblFinder:::.xgbtrain(X, y, nrounds = NULL, max_depth = 4, nthreads = 2)
                p <- predict(fit, X); sprintf("TRAIN_OK  AUC-proxy=%.3f", mean(p[y=="doublet"]) - mean(p[y=="real"])) },
              error = function(e) paste("TRAIN_ERROR:", conditionMessage(e)))
cat(r, "\n")
