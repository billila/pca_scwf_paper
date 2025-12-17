### create BE1 data ####

# download data from here: 
# https://figshare.com/articles/dataset/BE1_10XGenomics_count_matrices/23939481?file=42312711 


library(DropletUtils)
library(SingleCellExperiment)

data_dir <- "/home/ilaria/Downloads/BE1run12/BE1run12/"
dataset_list <- list.files(data_dir)
for (dataset in dataset_list) {
  current_data_dir <- file.path(data_dir, dataset)
  
  sce <- read10xCounts(
    samples = current_data_dir,
    sample.names = dataset,
    col.names = TRUE)
  
  new_name <- paste0("sce_", dataset)
  assign(new_name, sce)
  
}

sce <- cbind(sce_A549, `sce_CCL-185-IG`, sce_CRL5868, sce_DV90, 
             sce_HCC78, sce_HTB178, sce_PBMCs, sce_PC9)

# create RData object ####
save(sce, file = "BE1.RData")

# save in h5ad #####

counts(sce) <- as.matrix(counts(sce))

library(zellkonverter)
out_path <- tempfile(pattern = ".h5ad")

writeH5AD(sce, file = "BE1.h5ad")
