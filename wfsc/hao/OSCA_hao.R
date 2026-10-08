library(SingleCellExperiment)
library(scuttle)
library(scran)
library(scater)
library(mclust)
library(bluster)

time <- matrix(NA, 11, 1)
colnames(time) <- c("time_sec")
rownames(time) <- c("find_mit_gene", "filter", "normalization", "hvg",
                    "scaling", "PCA", "t-sne", "knn", "umap", "louvain", "leiden")

### data ####

sce <- readRDS("hao21_pbmc.rds")

#### 1. find mithocondial genes  ####
start_time <- Sys.time()
library(EnsDb.Hsapiens.v75)
chr.loc <- mapIds(EnsDb.Hsapiens.v75, keys=rownames(sce),
                  keytype="SYMBOL", column="SEQNAME")

is.mito <- which(chr.loc=="MT")
table(is.mito)

library(scuttle)
df <- perCellQCMetrics(sce, subsets=list(Mito=is.mito))

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
time_elapsed
print(paste("Time Elapsed:", time_elapsed))
time[1,1] <- time_elapsed

# filter data  ####
start_time <- Sys.time()
sce
end_time <- Sys.time()
time_elapsed <- end_time - start_time
time_elapsed
print(paste("Time Elapsed:", time_elapsed))
time[2,1] <- time_elapsed

# normalization ####
start_time <- Sys.time()
sce <- logNormCounts(sce)

end_time <- Sys.time()
time_elapsed <- end_time - start_time
time_elapsed
print(paste("Time Elapsed:", time_elapsed))
time[3,1] <- time_elapsed

# Identification of highly variable features (feature selection) ####
start_time <- Sys.time()
dec.sce <- modelGeneVar(sce)
fit.sce <- metadata(dec.sce)

hvg.sce.var <- getTopHVGs(dec.sce, n=1000)

end_time <- Sys.time()
time_elapsed <- end_time - start_time
print(paste("Time Elapsed:", time_elapsed))
time_elapsed
time[4,1] <- time_elapsed

save(list = "hvg.sce.var", file = "hao_bioc_hvg.RData")

# PCA ####
start_time <- Sys.time()
sce <- runPCA(sce, subset_row=hvg.sce.var)
sce

end_time <- Sys.time()
time_elapsed <- end_time - start_time
print(paste("Time Elapsed:", time_elapsed))
time_elapsed
time[6,1] <- time_elapsed

# t-sne ####
start_time <- Sys.time()
set.seed(100000)
sce <- runTSNE(sce, dimred="PCA")

end_time <- Sys.time()
time_elapsed <- end_time - start_time
print(paste("Time Elapsed:", time_elapsed))
time[7,1] <- time_elapsed


# knn ####

start_time <- Sys.time()
g <- buildSNNGraph(sce, use.dimred="PCA", k=50)
end_time <- Sys.time()
print(paste("Time Elapsed:", time_elapsed))
time_elapsed
time[8,1] <- time_elapsed

# umap ####
start_time <- Sys.time()
set.seed(1000000)
sce <- runUMAP(sce, dimred="PCA")

end_time <- Sys.time()
time_elapsed <- end_time - start_time
print(paste("Time Elapsed:", time_elapsed))
time_elapsed
time[9,1] <- time_elapsed


# louvain  ####
start_time <- Sys.time()
colLabels(sce) <- clusterCells(sce, use.dimred = "PCA",
                               BLUSPARAM = NNGraphParam(k = 50, cluster.fun = "louvain"))
end_time <- Sys.time()
time_elapsed <- end_time - start_time
time_elapsed
print(paste("Time Elapsed:", time_elapsed))

time[10,1] <- time_elapsed

table(colLabels(sce))
ARI <- adjustedRandIndex((sce$celltype.l1), colLabels(sce))
cat("Louvain Adjusted Rand Index:", ARI, "\n")



# leiden ####
start_time <- Sys.time()
colLabels(sce) <- clusterCells(sce, use.dimred = "PCA",
                               BLUSPARAM = NNGraphParam(k = 50, cluster.fun = "leiden"))
end_time <- Sys.time()
time_elapsed <- end_time - start_time
time_elapsed
print(paste("Time Elapsed:", time_elapsed))
time[11,1] <- time_elapsed

table(colLabels(sce))
ARI <- adjustedRandIndex((sce$celltype.l1), colLabels(sce))
ARI
cat("Leiden Adjusted Rand Index:", ARI, "\n")

time



save.image("/mnt/spca/pipeline_sc_2026/CITE_seq_161K/bioc/cite_bioc.RData")
