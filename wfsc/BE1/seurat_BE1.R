library(dplyr)
library(Seurat)
library(patchwork) 
library(SingleCellExperiment)

# save time usage #### 
time <- matrix(NA, 10, 1)
colnames(time) <- c("time_sec")
rownames(time) <- c("find_mit_gene", "filter", "normalization", "hvg", 
                    "scaling", "PCA", "t-sne", "umap", "louvain", "leiden")





load("BE1.RData")

# data ####
rownames(sce) <- rowData(sce)$Symbol
data <- as.Seurat(sce, counts = "counts", data = NULL, assay = NULL)
data

DefaultAssay(data) <- "originalexp"
DefaultAssay(data)

# find mitocondrial genes ####
start_time <- Sys.time()
data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^MT")
end_time <- Sys.time()
time_elapsed <- end_time - start_time
time_elapsed
print(paste("Time Elapsed:", time_elapsed))
time[1,1] <- time_elapsed

# filter data ####
start_time <- Sys.time()
VlnPlot(data, features = c("nFeature_originalexp", "nCount_originalexp", "percent.mt"), ncol = 3)
data <- subset(data, subset = nFeature_originalexp > 200 & nFeature_originalexp < 5000  & percent.mt < 5 &
                 nCount_originalexp < 25000)
end_time <- Sys.time()
time_elapsed <- end_time - start_time
time_elapsed
print(paste("Time Elapsed:", time_elapsed))
time[2,1] <- time_elapsed


# normalization ####
start_time <- Sys.time()
data <- NormalizeData(data, normalization.method = "LogNormalize")
end_time <- Sys.time()
time_elapsed <- end_time - start_time
time_elapsed
print(paste("Time Elapsed:", time_elapsed))
time[3,1] <- as.numeric(time_elapsed)

# Identification of highly variable features (feature selection) ####
start_time <- Sys.time()
data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 1000)
top10 <- head(VariableFeatures(data), 10)
top10
end_time <- Sys.time()
time_elapsed <- end_time - start_time
time_elapsed
print(paste("Time Elapsed:", time_elapsed))
time[4,1] <- as.numeric(time_elapsed)

hvg <- VariableFeatures(data)
#save(list = "hvg", file = "/mnt/spca/pipeline_sc/BE1/seurat/be1_seurat_hvg_symbol.RData")


# Scaling the data ####
start_time <- Sys.time()
all.genes <- rownames(data)
data <- ScaleData(data, features = all.genes)
end_time <- Sys.time()
time_elapsed <- end_time - start_time
time_elapsed
print(paste("Time Elapsed:", time_elapsed))
time[5,1] <- time_elapsed

# PCA ####
start_time <- Sys.time()
data <- RunPCA(data, features = VariableFeatures(object = data), verbose = FALSE)
end_time <- Sys.time()
time_elapsed <- end_time - start_time
time_elapsed
print(paste("Time Elapsed:", time_elapsed))
time[6,1] <- time_elapsed

# t-sne ####
start_time <- Sys.time()
data <- RunTSNE(data, 
                reduction = "pca", perplexity = 18)
end_time <- Sys.time()
time_elapsed <- end_time - start_time
time_elapsed
print(paste("Time Elapsed:", time_elapsed))
time[7,1] <- time_elapsed

# UMAP ####
start_time <- Sys.time()
data <- RunUMAP(data, dims = 1:50)
end_time <- Sys.time()
time_elapsed <- end_time - start_time
time_elapsed
print(paste("Time Elapsed:", time_elapsed))
time[8,1] <- time_elapsed


# louvain ####
start_time <- Sys.time()
data <- FindNeighbors(data, dims = 1:50, verbose = T)
data <- FindClusters(data, resolution = 0.2, algorithm = 1)
end_time <- Sys.time()
time_elapsed <- end_time - start_time
time_elapsed
print(paste("Time Elapsed:", time_elapsed))
time[9,1] <- time_elapsed

# cluster concordance louvain #####

num_clusters <- length(unique(data$seurat_clusters))
cat("Number of clusters:", num_clusters, "\n")
(table(data$seurat_clusters))
(table(data$Sample))


library(mclust)
ARI <- adjustedRandIndex((data$Sample), data$seurat_clusters)
cat("Louvain Adjusted Rand Index:", ARI, "\n")

# leiden ####
start_time <- Sys.time()
data <- FindNeighbors(data, dims = 1:50, verbose = T)
data <- FindClusters(data,  algorithm = 1, resolution = 0.02) # 4 leiden
end_time <- Sys.time()
time_elapsed <- end_time - start_time
time_elapsed
print(paste("Time Elapsed:", time_elapsed))
time[10,1] <- time_elapsed

# cluster concordance leiden #####

num_clusters <- length(unique(data$seurat_clusters))
cat("Number of clusters:", num_clusters, "\n")
(table(data$seurat_clusters))
(table(data$Sample))


library(mclust)
ARI <- adjustedRandIndex((data$Sample), data$seurat_clusters)
cat("Leiden Adjusted Rand Index:", ARI, "\n")

DimPlot(data, label = TRUE)

time
print(time)

# save.image("/mnt/spca/pipeline_sc/BE1/seurat/BE1_seurat.RData")