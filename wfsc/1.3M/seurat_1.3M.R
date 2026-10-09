library(dplyr)
library(Seurat)
library(patchwork) 
library(BPCells)
library(Seurat)
library(SeuratObject)
library(Azimuth)

options(future.globals.maxSize = 1e11)

# save time usage #### 
time <- matrix(NA, 10, 1)
colnames(time) <- c("time_sec")
rownames(time) <- c("find_mit_gene", "filter", "normalization", "hvg", 
                    "scaling", "PCA", "t-sne", "umap", "louvain", "leiden")

### pulizia db###


library(MultiAssayExperiment)
library(SingleCellMultiModal)
library(SingleCellExperiment)
library(ExperimentHub)

print("here")

brain.mat <- open_matrix_dir(dir = "/mnt/spca/pipeline_sc/brain_counts")
brain.mat <- Azimuth:::ConvertEnsembleToSymbol(mat = brain.mat, species = "mouse")

# Create Seurat Object
data <- CreateSeuratObject(counts = brain.mat)

print("brain.mat ok")


# find mitocondrial genes ####
start_time <- Sys.time()
data[["percent.mt"]] <- PercentageFeatureSet(data, pattern = "^MT-")
end_time <- Sys.time()
time_elapsed <- end_time - start_time
print(time_elapsed)
print(paste("Time Elapsed:", time_elapsed))
time[1,1] <- time_elapsed

# filter data sulle cellule, profondità di sequenziamento####
start_time <- Sys.time()
VlnPlot(data, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
data <- subset(data, subset = nFeature_RNA > 200 & nFeature_RNA < 2500  & percent.mt < 5 &
                 nCount_RNA < 4000)
end_time <- Sys.time()
time_elapsed <- end_time - start_time
print(paste("Time Elapsed:", time_elapsed))
time[2,1] <- time_elapsed
print(time_elapsed)


# normalization ####
start_time <- Sys.time()
data <- NormalizeData(data, normalization.method = "LogNormalize")
end_time <- Sys.time()
time_elapsed <- end_time - start_time
print(paste("Time Elapsed:", time_elapsed))
time[3,1] <- time_elapsed
print(time_elapsed)

# Identification of highly variable features (feature selection) ####
start_time <- Sys.time()
data <- FindVariableFeatures(data, selection.method = "vst", nfeatures = 1000)
top10 <- head(VariableFeatures(data), 10)
top10
end_time <- Sys.time()
time_elapsed <- end_time - start_time
print(paste("Time Elapsed:", time_elapsed))
time[4,1] <- time_elapsed
print(time_elapsed)

# Scaling the data ####
start_time <- Sys.time()
all.genes <- rownames(data)
data <- ScaleData(data, features = all.genes)
end_time <- Sys.time()
time_elapsed <- end_time - start_time
print(paste("Time Elapsed:", time_elapsed))
time[5,1] <- time_elapsed
print(time_elapsed)

# PCA ####
start_time <- Sys.time()
data <- RunPCA(data, features = VariableFeatures(object = data), verbose = FALSE)
end_time <- Sys.time()
time_elapsed <- end_time - start_time
print(paste("Time Elapsed:", time_elapsed))
time[6,1] <- time_elapsed
print(time_elapsed)

# t-sne ####
start_time <- Sys.time()
data <- RunTSNE(data, 
                reduction = "pca", perplexity = 18)
end_time <- Sys.time()
time_elapsed <- end_time - start_time
print(paste("Time Elapsed:", time_elapsed))
time[7,1] <- time_elapsed
print(time_elapsed)

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
data <- FindNeighbors(data, dims = 1:50, verbose = T)
data <- FindClusters(data, resolution = 0.2, algorithm = 1)
end_time <- Sys.time()
time_elapsed <- end_time - start_time
print(paste("Time Elapsed:", time_elapsed))
time[10,1] <- time_elapsed
print(time_elapsed)


# leiden ####
start_time <- Sys.time()
data <- FindNeighbors(data, dims = 1:50, verbose = T)
data <- FindClusters(data,  algorithm = 4, resolution = 0.2) # 4 leiden
end_time <- Sys.time()
time_elapsed <- end_time - start_time
print(paste("Time Elapsed:", time_elapsed))
time[11,1] <- time_elapsed
print(time_elapsed)

time 
print(time)


