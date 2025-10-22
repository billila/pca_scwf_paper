library(SingleCellExperiment)
library(zellkonverter)
library(BiocSingular)
library(here)
library(HDF5Array)
library(mbkmeans)
library(ggplot2)
library(scater)
library(scran)
library(BiocParallel)
library(DelayedMatrixStats)
library(rhdf5)

# HDF5 - random

#Create HDF5 data
# loading "1M_neurons_data.h5ad" obtained from python script "scanpytenx.py"
setwd("/mnt/spca/run_spca_2025/M18/time")



#### 100k

mat_100k <- HDF5Array("/mnt/spca/run_spca_2025/hdf5_sparse/sparse_100k/assays.h5", "assay001", as.sparse = TRUE)

time.start <- proc.time()

now <- format(Sys.time(), "%b%d%H%M%OS3")
out_name <- paste0(input, "_100k", "_", "ila",".out")

Rprof(filename = here("output",paste0("M1_Random", out_name)), append = FALSE, memory.profiling = TRUE)

invisible(random_pca <- BiocSingular::runPCA(mat_100k, rank = 50,
                                             center = TRUE, scale = FALSE,
                                             BSPARAM = BiocSingular::RandomParam(),
                                             #BPPARAM = BiocParallel::MulticoreParam(1)
))

Rprof(NULL)
time.end <- proc.time()
time100k_random<- time.end - time.start
time100k_random


# elapsed time in minute
time100k_random[3]/60
head(random_pca$x[,1:2])


#### 500k

mat_500k <- HDF5Array("/mnt/spca/run_spca_2025/hdf5_sparse/sparse_500k/assays.h5", "assay001", as.sparse = TRUE)

time.start <- proc.time()

now <- format(Sys.time(), "%b%d%H%M%OS3")
out_name <- paste0(input, "_500k", "_", "ila",".out")

Rprof(filename = here("output",paste0("M1_Random", out_name)), append = FALSE, memory.profiling = TRUE)


invisible(random_pca <- BiocSingular::runPCA(mat_500k, rank = 50,
                                             center = TRUE, scale = FALSE,
                                             BSPARAM = BiocSingular::RandomParam(),
                                             #BPPARAM = BiocParallel::MulticoreParam(1)
))

Rprof(NULL)

time.end <- proc.time()
time500k_random <- time.end - time.start
time500k_random
# elapsed time in minute
time500k_random[3]/60
head(random_pca$x[,1:2])


#### 1M

mat_1M <- HDF5Array("/mnt/spca/run_spca_2025/hdf5_sparse/sparse_1M/assays.h5", "assay001", as.sparse = TRUE)

time.start <- proc.time()

now <- format(Sys.time(), "%b%d%H%M%OS3")
out_name <- paste0(input, "_1000k", "_", "ila",".out")

Rprof(filename = here("output",paste0("M1_Random", out_name)), append = FALSE, memory.profiling = TRUE)

invisible(random_pca <- BiocSingular::runPCA(mat_1M, rank = 50,
                                             center = TRUE, scale = FALSE,
                                             BSPARAM = BiocSingular::RandomParam(),
                                             #BPPARAM = BiocParallel::MulticoreParam(1)
))

Rprof(NULL)

time.end <- proc.time()
time1000k_random <- time.end - time.start
time1000k_random

# elapsed time in minute
time1000k_random[3]/60
head(random_pca$x[,1:2])




#### 1.3M

mat_13M <- HDF5Array("/mnt/spca/run_spca_2025/hdf5_sparse/sparse_1.3M/assays.h5", "assay001", as.sparse = TRUE)

time.start <- proc.time()

now <- format(Sys.time(), "%b%d%H%M%OS3")
out_name <- paste0(input, "_1.3M", "_", "ila",".out")

Rprof(filename = here("output",paste0("M1_Random", out_name)), append = FALSE, memory.profiling = TRUE)

invisible(random_pca <- BiocSingular::runPCA(mat_13M, rank = 50,
                                             center = TRUE, scale = FALSE,
                                             BSPARAM = BiocSingular::RandomParam(),
                                             #BPPARAM = BiocParallel::MulticoreParam(1)
))

Rprof(NULL)

time.end <- proc.time()
time1.3M_random <- time.end - time.start
time1.3M_random

# elapsed time in minute
time1.3M_random[3]/60
head(random_pca$x[,1:2])
