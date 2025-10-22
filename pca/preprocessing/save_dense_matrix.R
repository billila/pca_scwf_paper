setwd("~/")

library(TENxBrainData)
tenx <- TENxBrainData()

rownames(tenx)<- rowData(tenx)$Symbol

# load the data preprocessed in python to filter the matrix
dati <- "1M_neurons_data.h5ad"
dati <- readH5AD(dati, use_hdf5 = TRUE)

tenx <- tenx[rownames(dati),]
assay(tenx)

# normalization 
tenx <- logNormCounts(tenx)

mat <- as.matrix(assay(tenx))
ncells <- c(100, 500, 1000)
cellidx <- lapply(ncells, function(n) sample(colnames(mat), n * 1000))
sapply(cellidx, length)
mat100k <- mat[,cellidx[[1]]]
mat500k <- mat[,cellidx[[2]]]
mat1M <- mat[,cellidx[[3]]]