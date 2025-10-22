library(rARPACK)
library(RSpectra)


load("data_M5.RData")

#### 100k

time.start <- proc.time()

now <- format(Sys.time(), "%b%d%H%M%OS3")
out_name <- paste0(input, "_100k", "_", "ila",".out")

Rprof(filename = here("output",paste0("M15_arpack", out_name)), append = FALSE, memory.profiling = TRUE)


invisible(arpack_pca <- RSpectra::svds( 
  A = mat100k, 
  k = 50, 
  opts = list(center = TRUE, scale = FALSE)
))

arpack_pca$d <- arpack_pca$d * sqrt(nrow(arpack_pca$u) - 1)
arpack_pca$x <- arpack_pca$u %*% diag(arpack_pca$d)

Rprof(NULL)
time.end <- proc.time()
time100k_arpack<- time.end - time.start
time100k_arpack

# elapsed time in minute
time100k_arpack[3]/60
head(arpack_pca$x[,1:2])


#### 500k

time.start <- proc.time()

now <- format(Sys.time(), "%b%d%H%M%OS3")
out_name <- paste0(input, "_500k", "_", "ila",".out")

Rprof(filename = here("output",paste0("M15_arpack", out_name)), append = FALSE, memory.profiling = TRUE)

invisible(arpack_pca <- RSpectra::svds( 
  A = mat500k, 
  k = 50, 
  opts = list(center = TRUE, scale = FALSE)
))

arpack_pca$d <- arpack_pca$d * sqrt(nrow(arpack_pca$u) - 1)
arpack_pca$x <- arpack_pca$u %*% diag(arpack_pca$d)

Rprof(NULL)
time.end <- proc.time()
time500k_arpack <- time.end - time.start
time500k_arpack
# elapsed time in minute
time500k_arpack[3]/60
head(arpack_pca$x[,1:2])


#### 1M


time.start <- proc.time()

now <- format(Sys.time(), "%b%d%H%M%OS3")
out_name <- paste0(input,"_1000k", "_", "ila",".out")

Rprof(filename = here("output",paste0("M15_arpack", out_name)), append = FALSE, memory.profiling = TRUE)

invisible(arpack_pca <- RSpectra::svds( 
  A = mat1M, 
  k = 50, 
  opts = list(center = TRUE, scale = FALSE)
))

arpack_pca$d <- arpack_pca$d * sqrt(nrow(arpack_pca$u) - 1)
arpack_pca$x <- arpack_pca$u %*% diag(arpack_pca$d)


Rprof(NULL)
time.end <- proc.time()
time1000k_arpack <- time.end - time.start
time1000k_arpack

# elapsed time in minute
time1000k_arpack[3]/60
head(arpack_pca$x[,1:2])




#### 1.3M

time.start <- proc.time()

now <- format(Sys.time(), "%b%d%H%M%OS3")
out_name <- paste0(input,"_1.3M", "_", "ila",".out")

Rprof(filename = here("output",paste0("M15_arpack", out_name)), append = FALSE, memory.profiling = TRUE)


invisible(arpack_pca <- RSpectra::svds( 
  A = mat, 
  k = 50, 
  opts = list(center = TRUE, scale = FALSE)
))

arpack_pca$d <- arpack_pca$d * sqrt(nrow(arpack_pca$u) - 1)
arpack_pca$x <- arpack_pca$u %*% diag(arpack_pca$d)

Rprof(NULL)
time.end <- proc.time()
time1.3M_arpack <- time.end - time.start
time1.3M_arpack

# elapsed time in minute
time1.3M_arpack[3]/60
head(arpack_pca$x[,1:2])





