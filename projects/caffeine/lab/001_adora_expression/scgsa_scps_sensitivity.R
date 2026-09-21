#!/usr/bin/env Rscript
# scPS stability check on a small gene set.
#
# The authors' scPS_2024.R calls RunPCA(..., npcs = 10, weight.by.var = FALSE)
# with `# approx = FALSE` left commented out. npcs = 10 cannot run on an
# 8-gene panel, and the truncated SVD (irlba) that Seurat uses by default warns
# "did not converge--results might be invalid" when asked for nearly all
# singular values of such a small matrix.
#
# This compares, for the minimal panel:
#   A) approx = TRUE  (Seurat default, irlba)  - what the main run used
#   B) approx = FALSE (exact svd)              - numerically safe here
# across several seeds, to see whether the scPS ranking is stable at all.

suppressPackageStartupMessages({ library(Seurat); library(Matrix) })

TMP <- Sys.getenv("SCGSA_WORK", unset = "scgsa_work")
DATA <- file.path(TMP, "data"); OUT <- file.path(TMP, "results")
atlas <- "hbca"

PANEL_MIN <- c("ADORA1", "GNAO1", "GNAI1", "KCNJ3", "KCNJ6",
               "CACNA1A", "CACNA1B", "ADCY5")

counts <- Read10X_h5(file.path(DATA, paste0(atlas, "_counts.h5")))
meta <- read.csv(file.path(DATA, paste0(atlas, "_cells.csv")),
                 stringsAsFactors = FALSE)
colnames(counts) <- meta$barcode
obj <- CreateSeuratObject(counts = counts)
obj <- NormalizeData(obj, verbose = FALSE)
genes <- intersect(PANEL_MIN, rownames(obj))
obj <- ScaleData(obj, features = genes, verbose = FALSE)

scps_once <- function(obj, genes, approx, seed) {
  set.seed(seed)
  npcs <- min(10, length(genes) - 1)
  o <- RunPCA(obj, features = genes, npcs = npcs, weight.by.var = FALSE,
              approx = approx, seed.use = seed, verbose = FALSE)
  emb <- Embeddings(o, "pca"); sdev <- Stdev(o, "pca")
  varExp <- (sdev^2) / sum(sdev^2)
  cumVar <- cumsum(sdev^2) / sum(sdev^2)
  idx <- which(cumVar > 0.5)
  maxCombPC <- if (length(idx) == 0) 1 else max(min(idx), 1)
  PC.X <- emb - min(emb)
  w <- sapply(seq_len(maxCombPC), function(i) PC.X[, i] * varExp[i])
  agg <- if (is.matrix(w) && ncol(w) > 1) sqrt(rowSums(w)) else sqrt(as.vector(w))
  dat <- GetAssayData(o, layer = "data")
  meanExpr <- colMeans(dat[rownames(dat) %in% genes, , drop = FALSE])
  tapply(agg * meanExpr[names(agg)], meta$group, mean)
}

rows <- list()
for (approx in c(TRUE, FALSE)) {
  for (seed in c(42, 1, 7)) {
    v <- tryCatch(scps_once(obj, genes, approx, seed),
                  error = function(e) { cat("FAIL approx=", approx, " seed=",
                                            seed, ": ", conditionMessage(e),
                                            "\n", sep = ""); NULL })
    if (is.null(v)) next
    rows[[length(rows) + 1]] <- data.frame(
      approx = approx, seed = seed, group = names(v),
      score = as.numeric(v), stringsAsFactors = FALSE)
    top <- names(sort(v, decreasing = TRUE))[1:3]
    cat(sprintf("approx=%-5s seed=%-3d top3: %s\n", approx, seed,
                paste(top, collapse = " | ")))
  }
}
out <- do.call(rbind, rows)
write.csv(out, file.path(OUT, "scps_sensitivity_hbca.csv"), row.names = FALSE)

# rank stability across runs
library(stats)
w <- reshape(out, idvar = "group", timevar = c("approx"), drop = "seed",
             direction = "wide")
cat("\nwrote", file.path(OUT, "scps_sensitivity_hbca.csv"), "\n")
