library(MultiAssayExperiment)
library(SingleCellMultiModal)
library(SingleCellExperiment)
library(scuttle)
library(AnnotationDbi)
library(scran)
library(scater)
library(mclust)
library(bluster)

# save time usage #### 
time <- matrix(NA, 10, 1)
colnames(time) <- c("time_sec")
rownames(time) <- c("find_mit_gene", "filter", "normalization", "hvg", 
                    "scaling", "PCA", "t-sne", "umap", "louvain", "leiden")

CITEseq(DataType="cord_blood", modes="*", dry.run=TRUE, version="1.0.0")

### data #### 

sce <- CITEseq(DataType="cord_blood", modes="*", dry.run=FALSE, version="1.0.0",
               DataClass="SingleCellExperiment")

colnames(sce)
rownames(sce) <- sub("^MOUSE_", "", rownames(sce))

gene_all <- rownames(sce)
gene_m <- gene_all[grep("HUMAN_", gene_all)]

sce_m <- sce[gene_m, ]
rownames(sce_m) <- sub("^HUMAN_", "", rownames(sce_m))
sce <- sce_m 


#### 1. find mithocondial genes  ####
start_time <- Sys.time()
library(EnsDb.Mmusculus.v79)
chr.loc <- mapIds(EnsDb.Mmusculus.v79, keys=rownames(sce),
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
print(paste("Time Elapsed:", time_elapsed))
time[1,1] <- time_elapsed

# filter data  ####
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
dec.sce <- modelGeneVar(sce)
fit.sce <- metadata(dec.sce)

hvg.sce.var <- getTopHVGs(dec.sce, n=1000)

end_time <- Sys.time()
time_elapsed <- end_time - start_time
print(paste("Time Elapsed:", time_elapsed))
time[4,1] <- time_elapsed

save(list = "hvg.sce.var", file = "/mnt/spca/pipeline_sc/cb/bioc/cb_bioc_hvg.RData")

# PCA ####
start_time <- Sys.time()
sce <- runPCA(sce, subset_row=hvg.sce.var)
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

table(colLabels(sce))
ARI <- adjustedRandIndex((sce$celltype), colLabels(sce))
cat("Louvain Adjusted Rand Index:", ARI, "\n")



# leiden ####
start_time <- Sys.time()
colLabels(sce) <- clusterCells(sce, use.dimred = "PCA",
                               BLUSPARAM = NNGraphParam(k = 50, cluster.fun = "leiden"))
end_time <- Sys.time()
time_elapsed <- end_time - start_time
print(paste("Time Elapsed:", time_elapsed))
time[10,1] <- time_elapsed

table(colLabels(sce))
ARI <- adjustedRandIndex((sce$celltype), colLabels(sce))
ARI
cat("Leiden Adjusted Rand Index:", ARI, "\n")

time

# save.image("cb_bioc.RData")
