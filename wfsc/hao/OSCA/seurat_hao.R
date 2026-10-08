### hao Seurat workflow ###

library(dplyr)
library(Seurat)
library(patchwork)
library(SingleCellExperiment)
library(mclust)


time <- matrix(NA, 11, 1)
colnames(time) <- c("time_sec")
rownames(time) <- c("find_mit_gene", "filter", "normalization", "hvg",
                    "scaling", "PCA", "t-sne", "knn", "umap", "louvain", "leiden")

sce <- readRDS("hao21_pbmc.rds")


# data ####
data <- as.Seurat(sce, counts = "counts", data = NULL, assay = NULL)
DefaultAssay(data) <- "originalexp"




# find mitocondrial genes ####
start_time <- Sys.time()
data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^MT-")
end_time <- Sys.time()
time_elapsed <- end_time - start_time
print(paste("Time Elapsed:", time_elapsed))
time[1,1] <- time_elapsed

# filter ####
start_time <- Sys.time()
VlnPlot(data, features = c("nFeature_originalexp", "nCount_originalexp", "percent.mt"), ncol = 3)
data <- subset(data, subset = nFeature_originalexp > 200 & nFeature_originalexp < 3000  & percent.mt < 5 &
                 nCount_originalexp < 15000)
end_time <- Sys.time()
time_elapsed <- end_time - start_time
print(paste("Time Elapsed:", time_elapsed))
time[2,1] <- time_elapsed


# normalization ####
start_time <- Sys.time()
data <- NormalizeData(data, normalization.method = "LogNormalize")
end_time <- Sys.time()
time_elapsed <- end_time - start_time
print(paste("Time Elapsed:", time_elapsed))
time[3,1] <- time_elapsed

# Identification of highly variable features (feature selection) ####
start_time <- Sys.time()
data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 1000)
top10 <- head(VariableFeatures(data), 10)
top10
end_time <- Sys.time()
time_elapsed <- end_time - start_time
print(paste("Time Elapsed:", time_elapsed))
time[4,1] <- time_elapsed

hvg <- VariableFeatures(data)
save(list = "hvg", file = "hao_seurat_hvg.RData")

# Scaling the data ####
start_time <- Sys.time()
all.genes <- rownames(data)
data <- ScaleData(data, features = all.genes)
end_time <- Sys.time()
time_elapsed <- end_time - start_time
print(paste("Time Elapsed:", time_elapsed))
time[5,1] <- time_elapsed

# PCA ####
start_time <- Sys.time()
data <- RunPCA(data, features = VariableFeatures(object = data), verbose = FALSE)
end_time <- Sys.time()
time_elapsed <- end_time - start_time
print(paste("Time Elapsed:", time_elapsed))
time[6,1] <- time_elapsed

# t-sne ####
start_time <- Sys.time()
data <- RunTSNE(data,
                reduction = "pca", perplexity = 18)
end_time <- Sys.time()
time_elapsed <- end_time - start_time
print(paste("Time Elapsed:", time_elapsed))
time[7,1] <- time_elapsed

# knn ####
start_time <- Sys.time()
data <- FindNeighbors(data, dims=1:50)
end_time <- Sys.time()
time_elapsed <- end_time - start_time
print(paste("Time Elapsed:", time_elapsed))
time[8,1] <- time_elapsed

# UMAP ####
start_time <- Sys.time()
data <- RunUMAP(data, dims = 1:50)
end_time <- Sys.time()
time_elapsed <- end_time - start_time
print(paste("Time Elapsed:", time_elapsed))
time[9,1] <- time_elapsed


# louvain ####
start_time <- Sys.time()

data <- FindClusters(data, resolution = 0.04, algorithm = 4)
end_time <- Sys.time()
time_elapsed <- end_time - start_time
print(paste("Time Elapsed:", time_elapsed))
time[10,1] <- time_elapsed

# cluster concordance louvain #####

num_clusters <- length(unique(data$seurat_clusters))
cat("Number of clusters:", num_clusters, "\n")
(table(data$seurat_clusters))
(table(data$celltype.l1))


library(mclust)
ARI <- adjustedRandIndex((data$celltype.l1), data$seurat_clusters)
cat("Louvain Adjusted Rand Index:", ARI, "\n")

# leiden ####
start_time <- Sys.time()

data <- FindClusters(data,  algorithm = 1, resolution = 0.04) 
end_time <- Sys.time()
time_elapsed <- end_time - start_time
print(paste("Time Elapsed:", time_elapsed))
time[11,1] <- time_elapsed

# cluster concordance leiden #####

num_clusters <- length(unique(data$seurat_clusters))
cat("Number of clusters:", num_clusters, "\n")
(table(data$seurat_clusters))
(table(data$celltype.l1))


library(mclust)
ARI <- adjustedRandIndex((data$celltype.l1), data$seurat_clusters)
cat("Leiden Adjusted Rand Index:", ARI, "\n")


time


##########################




save.image("/mnt/spca/pipeline_sc_2026/CITE_seq_161K/seurat/cite_seurat.RData")
