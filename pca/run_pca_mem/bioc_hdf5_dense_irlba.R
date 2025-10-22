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

### 100k

time.start <- proc.time()
sce <- loadHDF5SummarizedExperiment(dir = here("metodo2/data/subset/TENxBrainDataSE", "TENxBrainData_100k"), prefix="")

now <- format(Sys.time(), "%b%d%H%M%OS3")
out_name <- paste0(input, "_100k", "_", "ila",".out")

Rprof(filename = here("output", paste0("M1_Irlba", out_name)), append = FALSE, memory.profiling = TRUE)

invisible(random_pca <- BiocSingular::runPCA(assay(sce), rank = 50,
                                             center = TRUE, scale = FALSE,
                                             BSPARAM = BiocSingular::IrlbaParam(),
                                             #BPPARAM = BiocParallel::MulticoreParam(1)
))

Rprof(NULL)

time.end <- proc.time()
time100k_irlba<- time.end - time.start
time100k_irlba

# elapsed time in minute
time100k_irlba[3]/60
head(random_pca$x[,1:2])


#### 500k

time.start <- proc.time()
sce <- loadHDF5SummarizedExperiment(dir = here("metodo2/data/subset/TENxBrainDataSE", "TENxBrainData_500k"), prefix="")

now <- format(Sys.time(), "%b%d%H%M%OS3")
out_name <- paste0(input, "_500k", "_", "ila",".out")

Rprof(filename = here("output", paste0("M1_Irlba", out_name)), append = FALSE, memory.profiling = TRUE)

invisible(random_pca <- BiocSingular::runPCA(assay(sce), rank = 50,
                                             center = TRUE, scale = FALSE,
                                             BSPARAM = BiocSingular::IrlbaParam(),
                                             #BPPARAM = BiocParallel::MulticoreParam(1)
))

Rprof(NULL)

time.end <- proc.time()
time500k_irlba <- time.end - time.start
time500k_irlba
# elapsed time in minute
time500k_irlba[3]/60
head(random_pca$x[,1:2])


#### 1M


time.start <- proc.time()
sce <- loadHDF5SummarizedExperiment(dir = here("metodo2/data/subset/TENxBrainDataSE", "TENxBrainData_1000k"), prefix="")

now <- format(Sys.time(), "%b%d%H%M%OS3")
out_name <- paste0(input, "_1000k", "_", "ila",".out")

Rprof(filename = here("output", paste0("M1_Irlba", out_name)), append = FALSE, memory.profiling = TRUE)

invisible(random_pca <- BiocSingular::runPCA(assay(sce), rank = 50,
                                             center = TRUE, scale = FALSE,
                                             BSPARAM = BiocSingular::IrlbaParam(),
                                             #BPPARAM = BiocParallel::MulticoreParam(1)
))

Rprof(NULL)

time.end <- proc.time()
time1000k_irlba <- time.end - time.start
time1000k_irlba

# elapsed time in minute
time1000k_irlba[3]/60
head(random_pca$x[,1:2])




#### 1.3M

time.start <- proc.time()
sce <- loadHDF5SummarizedExperiment(dir = here("metodo2/data/subset/TENxBrainDataSE", "TENxBrainData_1.3M"), prefix="")

now <- format(Sys.time(), "%b%d%H%M%OS3")
out_name <- paste0(input,"_1.3M", "_", "ila",".out")

Rprof(filename = here("output", paste0("M1_Irlba", out_name)), append = FALSE, memory.profiling = TRUE)

invisible(random_pca <- BiocSingular::runPCA(assay(sce), rank = 50,
                                             center = TRUE, scale = FALSE,
                                             BSPARAM = BiocSingular::IrlbaParam(),
                                             #BPPARAM = BiocParallel::MulticoreParam(1)
))

Rprof(NULL)
time.end <- proc.time()
time1.3M_irlba <- time.end - time.start
time1.3M_irlba

# elapsed time in minute
time1.3M_irlba[3]/60
head(random_pca$x[,1:2])
