library(SingleCellExperiment)
library(Seurat)

setwd("/../GSE164378_RAW/")

meta1 <- read.csv("GSE164378_sc.meta.data_3P.csv")

mat1 <- ReadMtx(mtx = "GSM5008737_RNA_3P-matrix.mtx.gz",
                cells = "GSM5008737_RNA_3P-barcodes.tsv.gz",
                features = "GSM5008737_RNA_3P-features.tsv.gz", feature.column = 1)

hao <- SingleCellExperiment(assays = list(counts = mat1))
meta <- meta1[order(match(meta1$X, colnames(hao1))), ]
colData(hao) <- DataFrame(meta1)
hao1$technology <- "10x3"

saveRDS(sce, file = "hao21_pbmc.RData")


# save h5ad object #####
library(zellkonverter)
out_path <- tempfile(pattern = ".h5ad")

writeH5AD(sce, file = "hao21_pbmc.h5ad")