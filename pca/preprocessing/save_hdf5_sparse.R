library(HDF5Array)

# import my sparse data 

load("/mnt/spca/run_spca_2025/M5/time/data_M5.RData")


# 100k 

sce_100k <- sce[, colnames(sce) %in% mat100k@Dimnames[[2]]]
mat100k <- mat100k[rownames(sce_100k), colnames(sce_100k)]
counts(sce_100k) <- mat100k

saveHDF5SummarizedExperiment(sce_100k, 
                             dir = "/mnt/spca/run_spca_2025/hdf5_sparse/sparse_100k", 
                             prefix="", replace=TRUE, 
                             #chunkdim=c(dim(counts(sub))[1],1), 
                             level=NULL, verbose=FALSE)


# 500k 

sce_500k <- sce[, colnames(sce) %in% mat500k@Dimnames[[2]]]
mat500k <- mat500k[rownames(sce_500k), colnames(sce_500k)]
counts(sce_500k) <- mat500k

saveHDF5SummarizedExperiment(sce_500k, 
                             dir = "/mnt/spca/run_spca_2025/hdf5_sparse/sparse_500k", 
                             prefix="", replace=TRUE, 
                             #chunkdim=c(dim(counts(sub))[1],1), 
                             level=NULL, verbose=FALSE)


# 1M 

sce_1M <- sce[, colnames(sce) %in% mat1M@Dimnames[[2]]]
mat1M <- mat1M[rownames(sce_1M), colnames(sce_1M)]
counts(sce_1M) <- mat1M

saveHDF5SummarizedExperiment(sce_1M, 
                             dir = "/mnt/spca/run_spca_2025/hdf5_sparse/sparse_1M", 
                             prefix="", replace=TRUE, 
                             #chunkdim=c(dim(counts(sub))[1],1), 
                             level=NULL, verbose=FALSE)


# 1.3M 

sce_13M <- sce[, colnames(sce) %in% mat13M@Dimnames[[2]]]
mat13M <- mat[rownames(sce_13M), colnames(sce_13M)]
counts(sce_13M) <- mat13M

saveHDF5SummarizedExperiment(sce_13M, 
                             dir = "/mnt/spca/run_spca_2025/hdf5_sparse/sparse_1.3M/", 
                             prefix="", replace=TRUE, 
                             #chunkdim=c(dim(counts(sub))[1],1), 
                             level=NULL, verbose=FALSE)




# load an Large singleCellExperiment object saved before with saveHDF5SummarizedExperiment

# 100k
sce_100k <- loadHDF5SummarizedExperiment(dir = here("metodo2/data/subset/TENxBrainDataSE", "TENxBrainData_100k"), prefix="")

# 500k
sce_500k <- loadHDF5SummarizedExperiment(dir = here("metodo2/data/subset/TENxBrainDataSE", "TENxBrainData_500k"), prefix="")

# 1M
sce_1M <- loadHDF5SummarizedExperiment(dir = here("metodo2/data/subset/TENxBrainDataSE", "TENxBrainData_1000k"), prefix="")

# 1.3M
sce_1.3M <- loadHDF5SummarizedExperiment(dir = here("metodo2/data/subset/TENxBrainDataSE", "TENxBrainData_1306.127k"), prefix="")