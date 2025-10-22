library(rARPACK)
library(RSpectra)
library(HDF5Array)

# 100k
mat_100k <- HDF5Array("/mnt/spca/run_spca_2025/hdf5_sparse/sparse_100k/assays.h5", "assay001", as.sparse = TRUE)

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

pca <- simple_PCA(mat_100k, k = 50)

arpack_pca <- pca
arpack_pca$d <- arpack_pca$d * sqrt(nrow(arpack_pca$u) - 1)
arpack_pca$x <- arpack_pca$u %*% diag(arpack_pca$d) 

time.end <- proc.time()
time100k_arpack<- time.end - time.start
time100k_arpack

# elapsed time in minute
time100k_arpack[3]/60
head(arpack_pca$x[,1:2])

# 500k

mat_500k <- HDF5Array("/mnt/spca/run_spca_2025/hdf5_sparse/sparse_500k/assays.h5", "assay001", as.sparse = TRUE)

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

pca <- simple_PCA(mat_500k, k = 50)

arpack_pca <- pca
arpack_pca$d <- arpack_pca$d * sqrt(nrow(arpack_pca$u) - 1)
arpack_pca$x <- arpack_pca$u %*% diag(arpack_pca$d) 


time.end <- proc.time()
time500k_arpack<- time.end - time.start
time500k_arpack

# elapsed time in minute
time500k_arpack[3]/60
head(arpack_pca$x[,1:2])

# 1M

mat_1M <- HDF5Array("/mnt/spca/run_spca_2025/hdf5_sparse/sparse_1M/assays.h5", "assay001", as.sparse = TRUE)

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

pca <- simple_PCA(mat_1M, k = 50)

arpack_pca <- pca
arpack_pca$d <- arpack_pca$d * sqrt(nrow(arpack_pca$u) - 1)
arpack_pca$x <- arpack_pca$u %*% diag(arpack_pca$d) 




time.end <- proc.time()
time1000k_arpack<- time.end - time.start
time1000k_arpack

# elapsed time in minute
time1000k_arpack[3]/60
head(arpack_pca$x[,1:2])

# 1.3M

mat_13M <- HDF5Array("/mnt/spca/run_spca_2025/hdf5_sparse/sparse_1.3M/assays.h5", "assay001", as.sparse = TRUE)


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

pca <- simple_PCA(mat_13M, k = 50)

arpack_pca <- pca
arpack_pca$d <- arpack_pca$d * sqrt(nrow(arpack_pca$u) - 1)
arpack_pca$x <- arpack_pca$u %*% diag(arpack_pca$d) 


time.end <- proc.time()
time1.3M_arpack<- time.end - time.start
time1.3M_arpack

# elapsed time in minute
time1.3M_arpack[3]/60
head(arpack_pca$x[,1:2])
