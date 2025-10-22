import numpy as np
import pandas as pd
import scanpy as sc
import sys
import time
import os
import numpy as np
from scipy import sparse
import anndata as ad

sc.settings.verbosity  = 3
sc.settings.set_figure_params(dpi = 70)
sc.logging.print_versions()
sc._settings.ScanpyConfig(n_jobs = 1)


# 100k
os.chdir("/mnt/spca/run_spca_2025/M1/time/metodo2/data/subset/TENxBrainDataSE/TENxBrainData_100k/")



start_time = time.time()


file = "assays.h5"
str = "assay001"
adata=sc.read_hdf(file, str)
sa = sparse.csr_matrix(adata.X) 
adata2 = ad.AnnData(sa)

sc.tl.pca(adata2, svd_solver = "randomized")
sys.argv = time.time() - start_time
print("--- %s seconds ---" % int(sys.argv))

# 500k 
os.chdir("/mnt/spca/run_spca_2025/M1/time/metodo2/data/subset/TENxBrainDataSE/TENxBrainData_500k/")

start_time = time.time()

file = "assays.h5"
str = "assay001"
adata=sc.read_hdf(file, str)
sa = sparse.csr_matrix(adata.X) 
adata2 = ad.AnnData(sa)

sc.tl.pca(adata2, svd_solver = "randomized")
sys.argv = time.time() - start_time
print("--- %s seconds ---" % int(sys.argv))


# 1M 
os.chdir("/mnt/spca/run_spca_2025/M1/time/metodo2/data/subset/TENxBrainDataSE/TENxBrainData_1000k/")

start_time = time.time()

file = "assays.h5"
str = "assay001"
adata=sc.read_hdf(file, str)
sa = sparse.csr_matrix(adata.X) 
adata2 = ad.AnnData(sa)

sc.tl.pca(adata2, svd_solver = "randomized")
sys.argv = time.time() - start_time
print("--- %s seconds ---" % int(sys.argv))

# 1.3M
os.chdir("/mnt/spca/run_spca_2025/M1/time/metodo2/data/subset/TENxBrainDataSE/TENxBrainData_1.3M/")

start_time = time.time()

file = "assays.h5"
str = "assay001"
adata=sc.read_hdf(file, str)
sa = sparse.csr_matrix(adata.X) 
adata2 = ad.AnnData(sa)

sc.tl.pca(adata2, svd_solver = "randomized")
sys.argv = time.time() - start_time
print("--- %s seconds ---" % int(sys.argv))
