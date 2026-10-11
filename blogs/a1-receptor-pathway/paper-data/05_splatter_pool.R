# SSD, as the authors' "Rscript SSD and GS simulation.rmd" builds it:
#   Normal_Bcell_params <- splatEstimate(as.matrix(Normal_Bcell@assays$RNA@counts))
#   splatSimulate(setParams(params, batchCells = 2 * cellcounts, seed = repl), method = "single")
# Cells in a single-group Splatter simulation are drawn independently from the
# same fitted model, so instead of one simulation per run we simulate one pool
# of 4,000 cells; the browser draws 2n (and, for S3/S4, a second 2n) from it.
suppressPackageStartupMessages({ library(Seurat); library(splatter); library(Matrix) })
dir <- "/mnt/bulk/dl_bio/gse164381"
Normal_Bcell <- readRDS(file.path(dir, "naive_b.rds"))
cnt <- as.matrix(LayerData(Normal_Bcell, layer = "counts"))
t0 <- Sys.time()
Normal_Bcell_params <- splatEstimate(cnt)
cat("splatEstimate took", round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "min\n")
saveRDS(Normal_Bcell_params, file.path(dir, "naive_b_splat_params.rds"))
print(Normal_Bcell_params)
sim <- splatSimulate(setParams(Normal_Bcell_params, update = list(batchCells = 4000, seed = 1)),
                     method = "single", verbose = FALSE)
m <- as(counts(sim), "CsparseMatrix")
# Splatter genes are synthetic (means drawn from the fitted model, not tied to
# any real gene), so they keep Splatter's own names: Gene1, Gene2, ...
writeMM(m, file.path(dir, "ssd_counts.mtx"))
writeLines(rownames(m), file.path(dir, "ssd_genes.txt"))
writeLines(paste0("SSD_", colnames(m)), file.path(dir, "ssd_cells.txt"))
z <- rowMeans(m == 0)
cat("SSD pool:", ncol(m), "cells,", nrow(m), "genes, nnz", length(m@x),
    " median genes/cell", median(colSums(m > 0)), " zero-frac < 0.9:", sum(z < 0.9), "\n")
