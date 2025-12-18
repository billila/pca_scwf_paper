### 1.3M input data ####

library(SingleCellExperiment)

# load 1.3M data in your R session ####

library(TENxBrainData)
sce <- TENxBrainData()

# save h5ad object #####

library(zellkonverter)
out_path <- tempfile(pattern = ".h5ad")

writeH5AD(tenx, file = "brain_13M.h5ad")