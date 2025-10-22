library(rARPACK)
library(RSpectra)


load("dati_M2.RData")


#### 100k ####

time.start <- proc.time()


invisible(arpack_pca <- RSpectra::svds( 
  A = mat100k, 
  k = 50, 
  opts = list(center = TRUE, scale = FALSE)
))

arpack_pca$d <- arpack_pca$d * sqrt(nrow(arpack_pca$u) - 1)
arpack_pca$x <- arpack_pca$u %*% diag(arpack_pca$d)

time.end <- proc.time()
time100k_arpack<- time.end - time.start
time100k_arpack

# elapsed time in minute
time100k_arpack[3]/60
head(arpack_pca$x[,1:2])

#### 500k ####

time.start <- proc.time()


invisible(arpack_pca <- RSpectra::svds( 
  A = mat500k, 
  k = 50, 
  opts = list(center = TRUE, scale = FALSE)
))

arpack_pca$d <- arpack_pca$d * sqrt(nrow(arpack_pca$u) - 1)
arpack_pca$x <- arpack_pca$u %*% diag(arpack_pca$d)

time.end <- proc.time()
time500k_arpack<- time.end - time.start
time500k_arpack

# elapsed time in minute
time500k_arpack[3]/60
head(arpack_pca$x[,1:2])

#### 1M ####

time.start <- proc.time()


invisible(arpack_pca <- RSpectra::svds( 
  A = mat1M, 
  k = 50, 
  opts = list(center = TRUE, scale = FALSE)
))

arpack_pca$d <- arpack_pca$d * sqrt(nrow(arpack_pca$u) - 1)
arpack_pca$x <- arpack_pca$u %*% diag(arpack_pca$d)

time.end <- proc.time()
time1000k_arpack<- time.end - time.start
time1000k_arpack

# elapsed time in minute
time1000k_arpack[3]/60
head(arpack_pca$x[,1:2])

#### 1.3M ####

time.start <- proc.time()


invisible(arpack_pca <- RSpectra::svds( 
  A = mat, 
  k = 50, 
  opts = list(center = TRUE, scale = FALSE)
))

arpack_pca$d <- arpack_pca$d * sqrt(nrow(arpack_pca$u) - 1)
arpack_pca$x <- arpack_pca$u %*% diag(arpack_pca$d)

time.end <- proc.time()
time1.3M_arpack<- time.end - time.start
time1.3M_arpack

# elapsed time in minute
time1.3M_arpack[3]/60
head(arpack_pca$x[,1:2])



