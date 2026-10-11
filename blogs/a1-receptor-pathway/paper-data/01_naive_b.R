# Rebuild the input of Wang & Thakar 2024's RWSD: "Naive Normal B cell".
# Their script (Thakar-Lab/scPS, "Rscript RWSD and GS simulation.rmd") loads a
# Seurat object of naive B cells from the GSE164381 healthy donors, annotated
# "as described in the original paper". That object and its labels are not
# published, so we approximate it: the 4 healthy-donor matrices, standard
# Seurat QC/clustering, then keep clusters whose marker profile is naive B
# (IGHD/TCL1A/FCER2/IL4R high, CD27/IGHG/IGHA low).
suppressPackageStartupMessages({ library(Seurat); library(Matrix) })
dir <- "/mnt/bulk/dl_bio/gse164381"
healthy <- c("GSM5008746_2010113866", "GSM5008747_2010113873",
             "GSM5008748_2010113875", "GSM5008749_2010113904")
objs <- lapply(healthy, function(s) {
  m <- Read10X_h5(file.path(dir, paste0(s, "_G.filtered_feature_bc_matrix.h5")))
  o <- CreateSeuratObject(m, project = s, min.cells = 3, min.features = 200)
  o$donor <- s
  o
})
o <- merge(objs[[1]], objs[-1], add.cell.ids = sub("_.*", "", healthy))
o <- JoinLayers(o)
o[["percent.mt"]] <- PercentageFeatureSet(o, pattern = "^MT-")
cat("cells before QC:", ncol(o), "\n")
o <- subset(o, percent.mt < 15 & nFeature_RNA < 4000)
cat("cells after QC:", ncol(o), "\n")
o <- NormalizeData(o, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
o <- FindVariableFeatures(o, nfeatures = 2000, verbose = FALSE)
o <- ScaleData(o, verbose = FALSE)
o <- RunPCA(o, npcs = 20, verbose = FALSE)
o <- FindNeighbors(o, dims = 1:15, verbose = FALSE)
o <- FindClusters(o, resolution = 0.5, verbose = FALSE)

naive <- c("IGHD", "TCL1A", "FCER2", "IL4R", "BACH2", "IGHM")
memory <- c("CD27", "IGHG1", "IGHA1", "AIM2", "TNFRSF13B", "CD80")
other <- c("CD3E", "NKG7", "LYZ", "CD14", "PPBP", "HBB", "JCHAIN", "MZB1")
avg <- AverageExpression(o, features = c(naive, memory, other, "MS4A1", "CD19"),
                         layer = "data")$RNA
avg <- as.matrix(avg)
z <- colMeans(avg[naive, , drop = FALSE]) - colMeans(avg[memory, , drop = FALSE])
tab <- data.frame(cluster = colnames(avg), n = as.integer(table(Idents(o))[sub("^g", "", colnames(avg))]),
                  naive_minus_memory = round(z, 2), MS4A1 = round(avg["MS4A1", ], 2),
                  IGHD = round(avg["IGHD", ], 2), CD27 = round(avg["CD27", ], 2),
                  IGHG1 = round(avg["IGHG1", ], 2), CD3E = round(avg["CD3E", ], 2),
                  JCHAIN = round(avg["JCHAIN", ], 2))
print(tab[order(-tab$naive_minus_memory), ], row.names = FALSE)
saveRDS(o, file.path(dir, "healthy_clustered.rds"))
