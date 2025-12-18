### create sc_mix data ####

# download the data from here:
# https://github.com/LuyiTian/sc_mixology/blob/master/data/sincell_with_class_5cl.RData

# 5 cell line 
load("/home/ilaria/Downloads/sincell_with_class_5cl.RData")

sce_sc_10x_5cl_qc

# create RData object ####

save(sce_sc_10x_5cl_qc, file = "sc_mix.RData")


# save h5ad object #####

library(zellkonverter)
out_path <- tempfile(pattern = ".h5ad")

writeH5AD(sce_sc_10x_5cl_qc, file = "sc_mixolgy_10x_5cl.h5ad")
