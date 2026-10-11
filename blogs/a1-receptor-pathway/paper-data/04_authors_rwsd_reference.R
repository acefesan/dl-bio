# Run the authors' RWSD simulation (Thakar-Lab/scPS, "Comparison analysis/
# Rscript RWSD and GS simulation.rmd", scenarios 1 and 3, replicate 1) on our
# naive-B object, and save every random draw plus the resulting matrix so
# scgsa-site/validate/rwsd_check.js can replay them in the browser simulator.
# Changes from the original: Seurat 5 accessors (LayerData instead of
# @assays$RNA@data), our input path, no gene renaming, outputs written out.
suppressPackageStartupMessages({ library(Seurat); library(Matrix); library(jsonlite) })
dir <- "/mnt/bulk/dl_bio/gse164381"
out <- file.path(dir, "authors_reference"); dir.create(out, showWarnings = FALSE)
Normal_Bcell <- readRDS(file.path(dir, "naive_b.rds"))
cellcounts <- as.integer(Sys.getenv("CELLS", "200"))
Threshold <- 0.75
Signal_Level <- 1.2
set.seed(42)
sim.mc <- replicate(11, {X <- sample(1:length(colnames(Normal_Bcell)), (2 * cellcounts), replace = F)})
repl <- 1
D <- function(o) LayerData(o, layer = "data")

## ---- Scenario 1 (verbatim logic) ----
B_cell_subset <- Normal_Bcell[, colnames(Normal_Bcell)[sim.mc[, repl]]]
set.seed(42)
condition_1 <- sample(colnames(B_cell_subset), cellcounts, replace = F)
condition_2 <- setdiff(colnames(B_cell_subset), condition_1)
cond <- ifelse(colnames(B_cell_subset) %in% condition_1, "Condition_1", "Condition_2")
rows_zeros <- rowMeans(D(B_cell_subset) == 0)
Genes <- rownames(B_cell_subset)[rows_zeros < 1]
Signal <- rownames(B_cell_subset)[rows_zeros < Threshold]
NonSignal <- setdiff(Genes, Signal)
dat <- as.matrix(D(B_cell_subset))
dat <- dat[rownames(dat) %in% Genes, ]
dat[Signal, which(cond %in% "Condition_2")] <- dat[Signal, which(cond %in% "Condition_2")] * Signal_Level
cat("S1: detected", length(Genes), " signal", length(Signal), " non-signal", length(NonSignal), "\n")
cat("S1: median % cells expressing signal genes:",
    round(100 * median(1 - rows_zeros[Signal]), 1), "(paper reports 44.5 for RWSD)\n")
writeMM(Matrix(dat, sparse = TRUE), file.path(out, "s1_dat.mtx"))
write_json(list(cells = colnames(dat), cond = cond, genes = rownames(dat), signal = Signal),
           file.path(out, "s1.json"))

## ---- Scenario 3 (verbatim logic, Signal_Level as written) ----
B_cell_subset <- Normal_Bcell[, colnames(Normal_Bcell)[sim.mc[, repl]]]
cond3 <- c(rep("Condition_1", cellcounts), rep("Condition_2", cellcounts))
rows_zeros <- rowMeans(D(B_cell_subset) == 0)
Genes <- rownames(B_cell_subset)[rows_zeros < 1]
Signal <- rownames(B_cell_subset)[rows_zeros < Threshold]
dat <- as.matrix(D(B_cell_subset))
dat <- dat[rownames(dat) %in% Genes, ]
dat[Signal, which(cond3 %in% "Condition_2")] <- dat[Signal, which(cond3 %in% "Condition_2")] * Signal_Level
B2 <- Normal_Bcell[, colnames(Normal_Bcell)[sim.mc[, repl + 1]]]
rows_zeros2 <- rowMeans(D(B2) == 0)
Add_gene <- sample(rownames(B2)[rows_zeros2 < Threshold], 250, replace = F)
dat.1 <- as.matrix(D(B2))
dat.1 <- dat.1[rownames(dat.1) %in% Add_gene, ]
dat.1[, which(cond3 %in% "Condition_1")] <- 0
rownames(dat.1) <- paste0(rownames(dat.1), "add")
dat <- rbind(dat, dat.1)
cat("S3: matrix", nrow(dat), "x", ncol(dat), " added", length(Add_gene), "\n")
writeMM(Matrix(dat, sparse = TRUE), file.path(out, "s3_dat.mtx"))
write_json(list(cells = colnames(B_cell_subset), cells2 = colnames(B2), genes = rownames(dat),
                signal = Signal, added = sub("add$", "", rownames(dat.1))),
           file.path(out, "s3.json"))
