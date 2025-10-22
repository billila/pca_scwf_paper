library(rARPACK)
library(RSpectra)


load("data_M5.RData")


# 100k
sce <- loadHDF5SummarizedExperiment("/mnt/spca/run_spca_2025/M1/time/metodo2/data/subset/TENxBrainDataSE/TENxBrainData_100k")

now <- format(Sys.time(), "%b%d%H%M%OS3")
out_name <- paste0(input, "_100k", "_", "ila",".out")

Rprof(filename = here("output",paste0("M19_arpack", out_name)), append = FALSE, memory.profiling = TRUE)

time.start <- proc.time()

simple_PCA <- function(mat, k=25)
{
  stopifnot(length(dim(mat)) == 2)
  row_means <- rowMeans(mat)
  Ax <- function(x, args)
    (as.numeric(mat %*% x) - row_means * sum(x))
  Atx <- function(x, args)
    (as.numeric(x %*% mat) - as.vector(row_means %*% x))
  RSpectra::svds(Ax, Atrans=Atx, k=k, dim=dim(mat))
}

pca <- simple_PCA(assay(sce), k = 50)
arpack_pca <- pca
arpack_pca$d <- arpack_pca$d * sqrt(nrow(arpack_pca$u) - 1)
arpack_pca$x <- arpack_pca$u %*% diag(arpack_pca$d) 

Rprof(NULL)

time.end <- proc.time()
time100k_arpack<- time.end - time.start
time100k_arpack

# elapsed time in minute
time100k_arpack[3]/60
head(arpack_pca$x[,1:2])

# 500k
sce <- loadHDF5SummarizedExperiment("/mnt/spca/run_spca_2025/M1/time/metodo2/data/subset/TENxBrainDataSE/TENxBrainData_500k")

now <- format(Sys.time(), "%b%d%H%M%OS3")
out_name <- paste0(input, "_500k", "_", "ila",".out")

Rprof(filename = here("output",paste0("M19_arpack", out_name)), append = FALSE, memory.profiling = TRUE)


time.start <- proc.time()

simple_PCA <- function(mat, k=25)
{
  stopifnot(length(dim(mat)) == 2)
  row_means <- rowMeans(mat)
  Ax <- function(x, args)
    (as.numeric(mat %*% x) - row_means * sum(x))
  Atx <- function(x, args)
    (as.numeric(x %*% mat) - as.vector(row_means %*% x))
  RSpectra::svds(Ax, Atrans=Atx, k=k, dim=dim(mat))
}

pca <- simple_PCA(assay(sce), k = 50)
arpack_pca <- pca
arpack_pca$d <- arpack_pca$d * sqrt(nrow(arpack_pca$u) - 1)
arpack_pca$x <- arpack_pca$u %*% diag(arpack_pca$d) 

Rprof(NULL)

time.end <- proc.time()
time500k_arpack<- time.end - time.start
time500k_arpack

# elapsed time in minute
time500k_arpack[3]/60
head(arpack_pca$x[,1:2])

# 1M
sce <- loadHDF5SummarizedExperiment("/mnt/spca/run_spca_2025/M1/time/metodo2/data/subset/TENxBrainDataSE/TENxBrainData_1000k/")

now <- format(Sys.time(), "%b%d%H%M%OS3")
out_name <- paste0(input, "_1000k", "_", "ila",".out")

Rprof(filename = here("output",paste0("M19_arpack", out_name)), append = FALSE, memory.profiling = TRUE)


time.start <- proc.time()

simple_PCA <- function(mat, k=25)
{
  stopifnot(length(dim(mat)) == 2)
  row_means <- rowMeans(mat)
  Ax <- function(x, args)
    (as.numeric(mat %*% x) - row_means * sum(x))
  Atx <- function(x, args)
    (as.numeric(x %*% mat) - as.vector(row_means %*% x))
  RSpectra::svds(Ax, Atrans=Atx, k=k, dim=dim(mat))
}

pca <- simple_PCA(assay(sce), k = 50)
arpack_pca <- pca
arpack_pca$d <- arpack_pca$d * sqrt(nrow(arpack_pca$u) - 1)
arpack_pca$x <- arpack_pca$u %*% diag(arpack_pca$d) 

Rprof(NULL)

time.end <- proc.time()
time1000k_arpack<- time.end - time.start
time1000k_arpack


# elapsed time in minute
time1000k_arpack[3]/60
head(arpack_pca$x[,1:2])

# 1.3M
sce <- loadHDF5SummarizedExperiment("/mnt/spca/run_spca_2025/M1/time/metodo2/data/subset/TENxBrainDataSE/TENxBrainData_1.3M/")

now <- format(Sys.time(), "%b%d%H%M%OS3")
out_name <- paste0(input, "_1.3M", "_", "ila",".out")

Rprof(filename = here("output",paste0("M19_arpack", out_name)), append = FALSE, memory.profiling = TRUE)


time.start <- proc.time()

simple_PCA <- function(mat, k=25)
{
  stopifnot(length(dim(mat)) == 2)
  row_means <- rowMeans(mat)
  Ax <- function(x, args)
    (as.numeric(mat %*% x) - row_means * sum(x))
  Atx <- function(x, args)
    (as.numeric(x %*% mat) - as.vector(row_means %*% x))
  RSpectra::svds(Ax, Atrans=Atx, k=k, dim=dim(mat))
}

pca <- simple_PCA(assay(sce), k = 50)
arpack_pca <- pca
arpack_pca$d <- arpack_pca$d * sqrt(nrow(arpack_pca$u) - 1)
arpack_pca$x <- arpack_pca$u %*% diag(arpack_pca$d) 

Rprof(NULL)

time.end <- proc.time()
time1.3M_arpack<- time.end - time.start
time1.3M_arpack

# elapsed time in minute
time1.3M_arpack[3]/60
head(arpack_pca$x[,1:2])
