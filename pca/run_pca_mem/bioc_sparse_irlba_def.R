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

### 100k

time.start <- proc.time()

now <- format(Sys.time(), "%b%d%H%M%OS3")
out_name <- paste0(input, "_100k", "_", "ila",".out")

Rprof(filename = here("output",paste0("M5_Irlba", out_name)), append = FALSE, memory.profiling = TRUE)

invisible(random_pca <- BiocSingular::runPCA(mat100k, rank = 50,
                                             center = TRUE, scale = FALSE,
                                             BSPARAM = BiocSingular::IrlbaParam(deferred = TRUE),
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

now <- format(Sys.time(), "%b%d%H%M%OS3")
out_name <- paste0(input, "_500k", "_", "ila",".out")

Rprof(filename = here("output",paste0("M5_Irlba", out_name)), append = FALSE, memory.profiling = TRUE)
invisible(random_pca <- BiocSingular::runPCA(mat500k, rank = 50,
                                             center = TRUE, scale = FALSE,
                                             BSPARAM = BiocSingular::IrlbaParam(deferred = TRUE),
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

now <- format(Sys.time(), "%b%d%H%M%OS3")
out_name <- paste0(input, "_1000k", "_", "ila",".out")

Rprof(filename = here("output",paste0("M5_Irlba", out_name)), append = FALSE, memory.profiling = TRUE)

invisible(random_pca <- BiocSingular::runPCA(mat1M, rank = 50,
                                             center = TRUE, scale = FALSE,
                                             BSPARAM = BiocSingular::IrlbaParam(deferred = TRUE),
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

now <- format(Sys.time(), "%b%d%H%M%OS3")
out_name <- paste0(input, "_1.3M", "_", "ila",".out")

Rprof(filename = here("output",paste0("M5_Irlba", out_name)), append = FALSE, memory.profiling = TRUE)

invisible(random_pca <- BiocSingular::runPCA(mat, rank = 50,
                                             center = TRUE, scale = FALSE,
                                             BSPARAM = BiocSingular::IrlbaParam(deferred = TRUE),
                                             #BPPARAM = BiocParallel::MulticoreParam(1)
))

Rprof(NULL)
time.end <- proc.time()
time1.3M_irlba <- time.end - time.start
time1.3M_irlba

# elapsed time in minute
time1.3M_irlba[3]/60
head(random_pca$x[,1:2])



data100_random<- data.frame(dataset = "TENxBrain_100k",
                            ncells = 100000,
                            method = "random",
                            elapsed_time = time100k_random[3])

data500_random<- data.frame(dataset = "TENxBrain_500k",
                            ncells = 500000,
                            method = "random",
                            elapsed_time = time500k_random[3])

data1000_random<- data.frame(dataset = "TENxBrain_1000k",
                             ncells = 1000000,
                             method = "random",
                             elapsed_time = time1000k_random[3])

data1300_random<- data.frame(dataset = "TENxBrain_1.3M",
                             ncells = 1300000,
                             method = "random",
                             elapsed_time = time1.3M_random[3])

data_random <- rbind(data100_random, data500_random, data1000_random, data1300_random)
data_random

data100_exact<- data.frame(dataset = "TENxBrain_100k",
                           ncells = 100000,
                           method = "exact",
                           elapsed_time = time100k_exact[3])

data500_exact<- data.frame(dataset = "TENxBrain_500k",
                           ncells = 500000,
                           method = "exact",
                           elapsed_time = time500k_exact[3])

data1000_exact<- data.frame(dataset = "TENxBrain_1000k",
                            ncells = 1000000,
                            method = "exact",
                            elapsed_time = time1000k_exact[3])

data1300_exact<- data.frame(dataset = "TENxBrain_1.3M",
                            ncells = 1300000,
                            method = "exact",
                            elapsed_time = time1.3M_exact[3])
