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


#### 100k

time.start <- proc.time()
sce <- loadHDF5SummarizedExperiment(dir = here("metodo2/data/subset/TENxBrainDataSE", "TENxBrainData_100k"), prefix="")

invisible(random_pca <- BiocSingular::runPCA(assay(sce), rank = 50,
                                             center = TRUE, scale = FALSE,
                                             BSPARAM = BiocSingular::ExactParam(),
                                             #BPPARAM = BiocParallel::MulticoreParam(1)
))

time.end <- proc.time()
time100k_exact<- time.end - time.start
time100k_exact

# elapsed time in minute
time100k_exact[3]/60
head(random_pca$x[,1:2])


#### 500k

time.start <- proc.time()
sce <- loadHDF5SummarizedExperiment(dir = here("metodo2/data/subset/TENxBrainDataSE", "TENxBrainData_500k"), prefix="")

invisible(random_pca <- BiocSingular::runPCA(assay(sce), rank = 50,
                                             center = TRUE, scale = FALSE,
                                             BSPARAM = BiocSingular::ExactParam(),
                                             #BPPARAM = BiocParallel::MulticoreParam(1)
))

time.end <- proc.time()
time500k_exact <- time.end - time.start
time500k_exact
# elapsed time in minute
time500k_exact[3]/60
head(random_pca$x[,1:2])


#### 1M


time.start <- proc.time()
sce <- loadHDF5SummarizedExperiment(dir = here("metodo2/data/subset/TENxBrainDataSE", "TENxBrainData_1000k"), prefix="")

invisible(random_pca <- BiocSingular::runPCA(assay(sce), rank = 50,
                                             center = TRUE, scale = FALSE,
                                             BSPARAM = BiocSingular::ExactParam(),
                                             #BPPARAM = BiocParallel::MulticoreParam(1)
))

time.end <- proc.time()
time1000k_exact <- time.end - time.start
time1000k_exact

# elapsed time in minute
time1000k_exact[3]/60
head(random_pca$x[,1:2])




#### 1.3M

time.start <- proc.time()
sce <- loadHDF5SummarizedExperiment(dir = here("metodo2/data/subset/TENxBrainDataSE", "TENxBrainData_1.3M"), prefix="")

invisible(random_pca <- BiocSingular::runPCA(assay(sce), rank = 50,
                                             center = TRUE, scale = FALSE,
                                             BSPARAM = BiocSingular::ExactParam(),
                                             #BPPARAM = BiocParallel::MulticoreParam(1)
))
time.end <- proc.time()
time1.3M_exact <- time.end - time.start
time1.3M_exact

# elapsed time in minute
time1.3M_exact[3]/60
head(random_pca$x[,1:2])
