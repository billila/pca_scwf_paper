library(MultiAssayExperiment)
library(SingleCellExperiment)
library(scuttle)
library(AnnotationDbi)
library(scran)
library(scater)
library(mclust)
library(bluster)
library(scrapper)
library(DelayedArray)

# save time usage ####

load("BE1.RData")
table(sce$Sample)
time <- matrix(NA, 10, 1)
colnames(time) <- c("time_sec")
rownames(time) <- c("find_mit_gene", "filter", "normalization", "hvg", 
                    "scaling", "PCA", "t-sne", "umap", "louvain", "leiden")

library(EnsDb.Hsapiens.v75)

sce <- sce[!duplicated(rowData(sce)$Symbol), ]
assay(sce) <- DelayedArray(assay(sce))


assay(sce)
rownames(sce) <- rowData(sce)$Symbol
#### 1. find mithocondial genes  ####

start_time <- Sys.time()

is.mito <- grepl("^MT-", rownames(sce))
table(is.mito)

rna.qc.metrics <- computeRnaQcMetrics(assay(sce), subsets=list(mt=is.mito), 
                                      num.threads=nthreads)

rna.qc.thresholds <- suggestRnaQcThresholds(rna.qc.metrics)
rna.qc.filter <- filterRnaQcMetrics(rna.qc.thresholds, rna.qc.metrics)

filtered <- sce[,rna.qc.filter,drop=FALSE]


end_time <- Sys.time()
time_elapsed <- end_time - start_time
time_elapsed
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

size.factors <- centerSizeFactors(rna.qc.metrics$sum[rna.qc.filter])
normalized <- normalizeCounts(assay(filtered), size.factors)

end_time <- Sys.time()
time_elapsed <- end_time - start_time
time_elapsed
print(paste("Time Elapsed:", time_elapsed))
time[3,1] <- time_elapsed

assay(filtered, "normalized") <- normalized

# Identification of highly variable features (feature selection) ####
start_time <- Sys.time()
gene.var <- modelGeneVariances(assay(filtered, "normalized"), num.threads=nthreads)
hvg.sce.var <- chooseHighlyVariableGenes(gene.var$statistics$residuals, top = 1000)


end_time <- Sys.time()
time_elapsed <- end_time - start_time
time_elapsed
print(paste("Time Elapsed:", time_elapsed))
time[4,1] <- time_elapsed

#save(list = "hvg.sce.var", file = "/mnt/spca/pipeline_sc/BE1/scrapper/be1_scrapper_hvg.RData")

# PCA ####
start_time <- Sys.time()

pca <- runPca((assay(filtered, "normalized")[hvg.sce.var,]), num.threads=nthreads, number = 50)

end_time <- Sys.time()
time_elapsed <- end_time - start_time
time_elapsed
print(paste("Time Elapsed:", time_elapsed))
time[6,1] <- time_elapsed

reducedDim(filtered, "PCA") <- t(pca$components)
dim(reducedDim(filtered, "PCA")[,1:2])

pca_var_pct_scrap <- pca$variance.explained

# t-sne ####
start_time <- Sys.time()

tsne.out <- runTsne(pca$components, num.threads=nthreads)

end_time <- Sys.time()
time_elapsed <- end_time - start_time
time_elapsed
print(paste("Time Elapsed:", time_elapsed))
time[7,1] <- time_elapsed

reducedDim(filtered, "TSNE") <- tsne.out

plotTSNE(filtered, color_by = "Sample")


# umap ####
start_time <- Sys.time()
set.seed(1000000)

umap.out <- runUmap(pca$components, num.threads=nthreads)

end_time <- Sys.time()
time_elapsed <- end_time - start_time
print(paste("Time Elapsed:", time_elapsed))
time[8,1] <- time_elapsed

reducedDim(filtered, "UMAP") <- tsne.out

plotUMAP(filtered, color_by = "Sample")


# louvain  ####
start_time <- Sys.time()

snn.graph <- buildSnnGraph(pca$components, num.threads=nthreads)
clust.out <- clusterGraph(snn.graph, method = c("multilevel"), multilevel.resolution = 0.18)

end_time <- Sys.time()
time_elapsed <- end_time - start_time
print(paste("Time Elapsed:", time_elapsed))
time[9,1] <- time_elapsed


length(clust.out$membership)
length(filtered$Sample)
ARI <- adjustedRandIndex((filtered$Sample), clust.out$membership)
ARI
cat("Leiden Adjusted Rand Index:", ARI, "\n")


# leiden ####
start_time <- Sys.time()
snn.graph <- buildSnnGraph(pca$components, num.threads=nthreads)
clust.out <- clusterGraph(snn.graph, method = c("leiden"), leiden.resolution = 0.18)
end_time <- Sys.time()
time_elapsed <- end_time - start_time
print(paste("Time Elapsed:", time_elapsed))
time[10,1] <- time_elapsed

length(clust.out$membership)
length(filtered$Sample)
ARI <- adjustedRandIndex((filtered$Sample), clust.out$membership)
ARI
cat("Leiden Adjusted Rand Index:", ARI, "\n")

time 

# save.image("BE1_scrapper.RData")