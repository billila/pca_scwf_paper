library(MultiAssayExperiment)
library(SingleCellMultiModal)
library(SingleCellExperiment)
library(scuttle)
library(AnnotationDbi)
library(scran)
library(scater)
library(mclust)
library(bluster)
library(scrapper)
library(DelayedArray)

sce <- CITEseq(DataType="cord_blood", modes="*", dry.run=FALSE, version="1.0.0",
               DataClass="SingleCellExperiment")
sce
table(sce$celltype)
rownames(sce) <- sub("^MOUSE_", "", rownames(sce))

gene_all <- rownames(sce)
gene_m <- gene_all[grep("HUMAN_", gene_all)]

sce_m <- sce[gene_m, ]
rownames(sce_m) <- sub("^HUMAN_", "", rownames(sce_m))
sce <- sce_m 

time <- matrix(NA, 10, 1)
colnames(time) <- c("time_sec")
rownames(time) <- c("find_mit_gene", "filter", "normalization", "hvg", 
                    "scaling", "PCA", "t-sne", "umap", "louvain", "leiden")

library(EnsDb.Mmusculus.v79)

assay(sce) <- DelayedArray(assay(sce))

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
print(paste("Time Elapsed:", time_elapsed))
time[1,1] <- time_elapsed

# filter data ####
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
print(paste("Time Elapsed:", time_elapsed))
time[3,1] <- time_elapsed

assay(filtered, "normalized") <- normalized

# Identification of highly variable features (feature selection) ####
start_time <- Sys.time()
gene.var <- modelGeneVariances(assay(filtered, "normalized"), num.threads=nthreads)
hvg.sce.var <- chooseHighlyVariableGenes(gene.var$statistics$residuals, top = 1000)

# hvg.sce.var <- rownames(filtered[hvg.sce.var,])

end_time <- Sys.time()
time_elapsed <- end_time - start_time
print(paste("Time Elapsed:", time_elapsed))
time[4,1] <- time_elapsed

#save(list = "hvg.sce.var", file = "/mnt/spca/pipeline_sc/cb/scrapper/cb_scrapper_hvg.RData")

# PCA ####
start_time <- Sys.time()

pca <- runPca((assay(filtered, "normalized")[hvg.sce.var,]), num.threads=nthreads, number = 50)

end_time <- Sys.time()
time_elapsed <- end_time - start_time
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
print(paste("Time Elapsed:", time_elapsed))
time[7,1] <- time_elapsed

reducedDim(filtered, "TSNE") <- tsne.out

plotTSNE(filtered, color_by = "celltype")

pca_var_pct_scrap <- pca$variance.explaine

# umap ####
start_time <- Sys.time()
set.seed(1000000)

umap.out <- runUmap(pca$components, num.threads=nthreads)

end_time <- Sys.time()
time_elapsed <- end_time - start_time
print(paste("Time Elapsed:", time_elapsed))
time[8,1] <- time_elapsed

reducedDim(filtered, "UMAP") <- umap.out

plotUMAP(filtered, color_by = "celltype")

# louvain  ####
start_time <- Sys.time()

snn.graph <- buildSnnGraph(pca$components, num.threads=nthreads)
clust.out <- clusterGraph(snn.graph, method = c("multilevel"), multilevel.resolution = 0.18)
end_time <- Sys.time()
time_elapsed <- end_time - start_time
print(paste("Time Elapsed:", time_elapsed))
time[9,1] <- time_elapsed


length(clust.out$membership)
length(filtered$celltype)
ARI <- adjustedRandIndex((filtered$celltype), clust.out$membership)
ARI
cat("Leiden Adjusted Rand Index:", ARI, "\n")


# leiden ####
start_time <- Sys.time()
snn.graph <- buildSnnGraph(pca$components, num.threads=nthreads)
clust.out <- clusterGraph(snn.graph, method = c("leiden"), leiden.resolution = 0.20)
end_time <- Sys.time()
time_elapsed <- end_time - start_time
print(paste("Time Elapsed:", time_elapsed))
time[10,1] <- time_elapsed

#table(colLabels(filtered))
length(clust.out$membership)
length(filtered$celltype)
ARI <- adjustedRandIndex((filtered$celltype), clust.out$membership)
ARI
cat("Leiden Adjusted Rand Index:", ARI, "\n")

time 

# save.image("cb_scrapper.RData")
