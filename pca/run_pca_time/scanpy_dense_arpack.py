import numpy as np
import pandas as pd
import scanpy as sc
import sys
import time
import os

sc.settings.verbosity  = 3
sc.settings.set_figure_params(dpi = 70)
sc.logging.print_versions()


os.chdir("/mnt/spca/run_spca_2025/M1/time/metodo2/data/subset/TENxBrainDataSE/TENxBrainData_100k/")


start_time = time.time()

file = "assays.h5"
str = "assay001"
adata=sc.read_hdf(file, str)


sc.tl.pca(adata, svd_solver = "arpack", use_highly_variable = False, n_comps = 50)
sys.argv = time.time() - start_time
print("--- %s seconds ---" % int(sys.argv))


os.chdir("/mnt/spca/run_spca_2025/M1/time/metodo2/data/subset/TENxBrainDataSE/TENxBrainData_500k/")


start_time = time.time()

file = "assays.h5"
str = "assay001"
adata=sc.read_hdf(file, str)

sc.tl.pca(adata, svd_solver = "arpack", use_highly_variable = False, n_comps = 50)
sys.argv = time.time() - start_time
print("--- %s seconds ---" % int(sys.argv))


os.chdir("/mnt/spca/run_spca_2025/M1/time/metodo2/data/subset/TENxBrainDataSE/TENxBrainData_1000k/")


start_time = time.time()

file = "assays.h5"
str = "assay001"
adata=sc.read_hdf(file, str)

sc.tl.pca(adata, svd_solver = "arpack", use_highly_variable = False, n_comps = 50)
sys.argv = time.time() - start_time
print("--- %s seconds ---" % int(sys.argv))


os.chdir("/mnt/spca/run_spca_2025/M1/time/metodo2/data/subset/TENxBrainDataSE/TENxBrainData_1.3M/")


start_time = time.time()

file = "assays.h5"
str = "assay001"
adata=sc.read_hdf(file, str)

sc.tl.pca(adata, svd_solver = "arpack", use_highly_variable = False, n_comps = 50)
sys.argv = time.time() - start_time
print("--- %s seconds ---" % int(sys.argv))
