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

load("data_M5.RData")

#### 100k

time.start <- proc.time()

now <- format(Sys.time(), "%b%d%H%M%OS3")
out_name <- paste0(input,"_100k", "_", "ila",".out")

Rprof(filename = here("output",paste0("M5_Exact", out_name)), append = FALSE, memory.profiling = TRUE)

invisible(random_pca <- BiocSingular::runPCA(mat100k, rank = 50,
                                             center = TRUE, scale = FALSE,
                                             BSPARAM = BiocSingular::ExactParam(deferred = TRUE),
                                             #BPPARAM = BiocParallel::MulticoreParam(1)
))

Rprof(NULL)

time.end <- proc.time()
time100k_exact<- time.end - time.start
time100k_exact

# elapsed time in minute
time100k_exact[3]/60
head(random_pca$x[,1:2])


#### 500k

time.start <- proc.time()

now <- format(Sys.time(), "%b%d%H%M%OS3")
out_name <- paste0(input,"_500k", "_", "ila",".out")

Rprof(filename = here("output",paste0("M5_Exact", out_name)), append = FALSE, memory.profiling = TRUE)

invisible(random_pca <- BiocSingular::runPCA(mat500k, rank = 50,
                                             center = TRUE, scale = FALSE,
                                             BSPARAM = BiocSingular::ExactParam(deferred = TRUE),
                                             #BPPARAM = BiocParallel::MulticoreParam(1)
))

Rprof(NULL)

time.end <- proc.time()
time500k_exact <- time.end - time.start
time500k_exact
# elapsed time in minute
time500k_exact[3]/60
head(random_pca$x[,1:2])


#### 1M


time.start <- proc.time()

now <- format(Sys.time(), "%b%d%H%M%OS3")
out_name <- paste0(input, "_1000k", "_", "ila",".out")

Rprof(filename = here("output",paste0("M5_Exact", out_name)), append = FALSE, memory.profiling = TRUE)

invisible(random_pca <- BiocSingular::runPCA(mat1M, rank = 50,
                                             center = TRUE, scale = FALSE,
                                             BSPARAM = BiocSingular::ExactParam(deferred = TRUE),
                                             #BPPARAM = BiocParallel::MulticoreParam(1)
))

Rprof(NULL)
time.end <- proc.time()
time1000k_exact <- time.end - time.start
time1000k_exact

# elapsed time in minute
time1000k_exact[3]/60
head(random_pca$x[,1:2])




#### 1.3M

time.start <- proc.time()

now <- format(Sys.time(), "%b%d%H%M%OS3")
out_name <- paste0(input, "_1.3M", "_", "ila",".out")

Rprof(filename = here("output",paste0("M5_Exact", out_name)), append = FALSE, memory.profiling = TRUE)

invisible(random_pca <- BiocSingular::runPCA(mat, rank = 50,
                                             center = TRUE, scale = FALSE,
                                             BSPARAM = BiocSingular::ExactParam(deferred = TRUE),
                                             #BPPARAM = BiocParallel::MulticoreParam(1)
))

Rprof(NULL)

time.end <- proc.time()
time1.3M_exact <- time.end - time.start
time1.3M_exact

# elapsed time in minute
time1.3M_exact[3]/60
head(random_pca$x[,1:2])
