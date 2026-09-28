#!/usr/bin/env Rscript
# Score the dataset exported by export.js with the real R packages, and
# compare against the JS engine's scores. Uses the same parameters irGSEA
# passes in the paper's comparison script.
# usage: Rscript validate.R <dir> <jasmine_functions.R>
suppressPackageStartupMessages({
  library(Matrix); library(Seurat); library(AUCell); library(UCell); library(GSVA)
})
args <- commandArgs(trailingOnly = TRUE)
dir <- if (length(args) >= 1) args[1] else "/tmp/scgsa-validate"
source(if (length(args) >= 2) args[2] else "jasmine_functions.R")

counts <- as.matrix(read.csv(file.path(dir, "counts.csv"), row.names = 1))
data_js <- as.matrix(read.csv(file.path(dir, "data_js.csv"), row.names = 1))
set <- readLines(file.path(dir, "set.txt"))
p <- jsonlite::fromJSON(file.path(dir, "params.json"))
js <- read.csv(file.path(dir, "scores_js.csv"), row.names = 1)

obj <- suppressWarnings(CreateSeuratObject(counts = as(counts, "dgCMatrix")))
# The simulator injects signal into the log-normalised matrix (as the paper's
# real-world simulation does), so every method is scored on that matrix.
obj <- SetAssayData(obj, layer = "data", new.data = as(data_js, "dgCMatrix"))
expr <- data_js
if (identical(p$inject, "counts")) {
  chk <- as.matrix(GetAssayData(NormalizeData(obj, verbose = FALSE), layer = "data"))
  cat(sprintf("normalisation max |R - JS| = %.2e\n", max(abs(chk - data_js))))
}

r <- list()

o <- suppressWarnings(AddModuleScore(obj, features = list(set), ctrl = p$ctrl, name = "AMS"))
r$AddModuleScore <- o$AMS1

rk <- AUCell_buildRankings(expr, plotStats = FALSE, verbose = FALSE)
r$AUCell <- getAUC(AUCell_calcAUC(list(s = set), rk, aucMaxRank = p$aucMaxRank,
                                  verbose = FALSE))[1, ]

uc <- ScoreSignatures_UCell(expr, features = list(s = set), maxRank = p$ucellMaxRank,
                            name = "", w_neg = 1)
r$UCell <- uc[, 1]

r$ssGSEA <- GSVA::gsva(GSVA::ssgseaParam(expr, list(s = set), normalize = FALSE),
                       verbose = FALSE)[1, ]

jas <- JASMINE(expr, set, method = "oddsratio")
r$JASMINE <- setNames(jas$JAS_Scores, jas$SampleID)

r$SCSE <- 100 * colSums(expr[set, , drop = FALSE]) / colSums(expr)

k <- length(set)
npcs <- min(10, k - 1)
o <- ScaleData(obj, features = set, verbose = FALSE)
# Exact SVD for stable numbers, but keep only the first npcs standard
# deviations: that is what the authors' default approx = TRUE (irlba) path
# stores, so "50% of variance" means 50% of the top-npcs variance.
o <- RunPCA(o, features = set, npcs = npcs, weight.by.var = FALSE, approx = FALSE,
            verbose = FALSE)
emb <- Embeddings(o, "pca"); sdev <- Stdev(o, "pca")[seq_len(npcs)]
ve <- sdev^2 / sum(sdev^2); cv <- cumsum(ve)
m <- min(which(cv > 0.5))
scps_from <- function(emb) {
  PCX <- emb - min(emb)
  agg <- sqrt(rowSums(PCX[, 1:m, drop = FALSE] %*% diag(ve[1:m], nrow = m)))
  agg * colMeans(expr[set, , drop = FALSE])
}
r$scPS <- scps_from(emb)
# Eigenvector signs are arbitrary. Re-express R's PCs under the JS convention
# (largest-magnitude loading positive) to isolate sign from everything else.
load <- Loadings(o, "pca")
flip <- apply(load, 2, function(v) if (v[which.max(abs(v))] < 0) -1 else 1)
r$scPS_signAligned <- scps_from(sweep(emb, 2, flip, `*`))
js$scPS_signAligned <- js$scPS
cat(sprintf("scPS R: npcs=%d maxCombPC=%d | JS: npcs=%d maxComb=%d\n",
            npcs, m, p$scps$npcs, p$scps$maxComb))

cat(sprintf("\n%-16s %10s %12s %14s\n", "method", "pearson", "spearman", "max|diff|"))
for (nm in names(r)) {
  a <- as.numeric(r[[nm]][rownames(js)]); b <- js[[nm]]
  cat(sprintf("%-16s %10.6f %12.6f %14.3e\n", nm, cor(a, b), cor(a, b, method = "spearman"),
              max(abs(a - b))))
}
