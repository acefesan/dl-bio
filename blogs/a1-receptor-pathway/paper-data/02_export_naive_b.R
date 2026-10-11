# Keep the naive-B clusters from 01_naive_b.R and save:
#  - naive_b.rds: Seurat object, LogNormalized (stand-in for the authors'
#    "Naive Normal B cell.RData"), used by the authors' simulation code.
#  - naive_b_counts.mtx/genes/cells: raw counts for the browser export.
# Clusters chosen from the marker table printed by 01 (IGHD/TCL1A high,
# CD27/IGHG1/CD3E/JCHAIN low): 0, 1, 3, 4, 7.
suppressPackageStartupMessages({ library(Seurat); library(Matrix) })
dir <- "/mnt/bulk/dl_bio/gse164381"
o <- readRDS(file.path(dir, "healthy_clustered.rds"))
Normal_Bcell <- subset(o, idents = c("0", "1", "3", "4", "7"))
Normal_Bcell <- NormalizeData(Normal_Bcell, normalization.method = "LogNormalize",
                              scale.factor = 10000, verbose = FALSE)
cat("naive B cells:", ncol(Normal_Bcell), " genes:", nrow(Normal_Bcell), "\n")
print(table(Normal_Bcell$donor))
saveRDS(Normal_Bcell, file.path(dir, "naive_b.rds"))
cnt <- LayerData(Normal_Bcell, layer = "counts")
writeMM(cnt, file.path(dir, "naive_b_counts.mtx"))
writeLines(rownames(cnt), file.path(dir, "naive_b_genes.txt"))
writeLines(colnames(cnt), file.path(dir, "naive_b_cells.txt"))
cat("nnz:", length(cnt@x), " max count:", max(cnt@x), "\n")
z <- rowMeans(LayerData(Normal_Bcell, layer = "data") == 0)
cat("genes detected:", sum(z < 1), " zero-frac < 0.75 (authors' 'Signal'):", sum(z < 0.75), "\n")
