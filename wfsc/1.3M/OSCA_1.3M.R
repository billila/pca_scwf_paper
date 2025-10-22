### cite_seq ###

### pulizia db###


library(MultiAssayExperiment)
library(SingleCellMultiModal)
library(SingleCellExperiment)
library(scuttle)
library(AnnotationDbi)
library(scran)
library(scater)
library(mclust)
library(bluster)
library(HDF5Array)
library(DelayedArray)
library(BiocSingular)

# save time usage #### 
time <- matrix(NA, 10, 1)
colnames(time) <- c("time_sec")
rownames(time) <- c("find_mit_gene", "filter", "normalization", "hvg", 
                    "scaling", "PCA", "t-sne", "umap", "louvain", "leiden")

library(TENxBrainData)
sce <- TENxBrainData()

#### 1. find mithocondial genes  ####
rownames(sce) <- rowData(sce)$Symbol
start_time <- Sys.time()
library(EnsDb.Mmusculus.v79)
chr.loc <- mapIds(EnsDb.Mmusculus.v79, keys=rownames(sce),
                  keytype="SYMBOL", column="SEQNAME")

is.mito <- which(chr.loc=="MT")
table(is.mito)

library(scuttle)
df <- perCellQCMetrics(sce, subsets=list(Mito=is.mito))

# include them in the object
colData(sce) <- cbind(colData(sce), df)
colData(sce)

table(df$sum < 10000)

table(df$subsets_Mito_percent > 5)

summary(df$detected)

summary(df$subsets_Mito_percent)

reasons <- perCellQCFilters(df, sub.fields="subsets_Mito_percent")
colSums(as.matrix(reasons))
summary(reasons$discard)
sce$discard <- reasons$discard

sce <- sce[,!sce$discard]
sce

end_time <- Sys.time()
time_elapsed <- end_time - start_time
print(paste("Time Elapsed:", time_elapsed))
time[1,1] <- time_elapsed

# filter data sulle cellule, profondità di sequenziamento ####
start_time <- Sys.time()
sce 
end_time <- Sys.time()
time_elapsed <- end_time - start_time
print(paste("Time Elapsed:", time_elapsed))
time[2,1] <- time_elapsed

# normalization ####
start_time <- Sys.time()
sce <- logNormCounts(sce)

end_time <- Sys.time()
time_elapsed <- end_time - start_time
print(paste("Time Elapsed:", time_elapsed))
time[3,1] <- time_elapsed

# Identification of highly variable features (feature selection) ####
start_time <- Sys.time()

top_hvg <- rowVars(counts(sce))

top_hvg <- top_hvg[order(top_hvg, decreasing = TRUE)]
gene_hvg <- names(top_hvg)[1:1000]
sce <- sce[ rownames(sce) %in% gene_hvg,]

end_time <- Sys.time()
time_elapsed <- end_time - start_time
print(paste("Time Elapsed:", time_elapsed))
time[4,1] <- time_elapsed

#save(list = "gene_hvg", file = "/mnt/spca/pipeline_sc/sc_mix/seurat/sc_mix_seurat_hvg.RData")

# PCA ####
start_time <- Sys.time()
sce <- BiocSingular::runPCA(sce, 50, BPPARAM=BiocParallel::MulticoreParam(1), BSPARAM=BiocSingular::RandomParam())
sce

end_time <- Sys.time()
time_elapsed <- end_time - start_time
print(paste("Time Elapsed:", time_elapsed))
time[6,1] <- time_elapsed

# t-sne ####
start_time <- Sys.time()
set.seed(100000)
sce <- runTSNE(sce, dimred="PCA")

end_time <- Sys.time()
time_elapsed <- end_time - start_time
print(paste("Time Elapsed:", time_elapsed))
time[7,1] <- time_elapsed

# umap ####
start_time <- Sys.time()
set.seed(1000000)
sce <- runUMAP(sce, dimred="PCA")

end_time <- Sys.time()
time_elapsed <- end_time - start_time
print(paste("Time Elapsed:", time_elapsed))
time[8,1] <- time_elapsed


# louvain  ####
start_time <- Sys.time()
colLabels(sce) <- clusterCells(sce, use.dimred = "PCA",
                               BLUSPARAM = NNGraphParam(k = 50, cluster.fun = "louvain"))
end_time <- Sys.time()
time_elapsed <- end_time - start_time
print(paste("Time Elapsed:", time_elapsed))
time[9,1] <- time_elapsed


# leiden ####
start_time <- Sys.time()
colLabels(sce) <- clusterCells(sce, use.dimred = "PCA",
                               BLUSPARAM = NNGraphParam(k = 50, cluster.fun = "leiden"))
end_time <- Sys.time()
time_elapsed <- end_time - start_time
print(paste("Time Elapsed:", time_elapsed))
time[10,1] <- time_elapsed

