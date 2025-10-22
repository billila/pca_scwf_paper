import cupy as cp
import rmm
import numpy as np
import scanpy as sc
import anndata
import time
import os
import scipy
import sys
import time
from cuml.decomposition import PCA
import pandas as pd
import GPUtil
n_components = 50


# 100k jacobi
adata = sc.read_csv("/media/volume/sdb/mat100k.csv", delimiter=',', first_column_names= True, dtype='float32')

GPUtil.showUtilization()
from pynvml import *
nvmlInit()
handle = nvmlDeviceGetHandleByIndex(0)

start_time = time.time()
adata.obsm["X_pca"] = PCA(n_components=n_components, svd_solver = "jacobi", output_type="numpy").fit_transform(adata.X)
p = PCA(n_components=n_components, output_type="numpy")

sys.argv = time.time() - start_time
print("--- %s seconds ---" % int(sys.argv))
GPUtil.showUtilization()
info = nvmlDeviceGetMemoryInfo(handle)
print("Total memory:", info.total)
print("Free memory:", info.free)
print("Used memory:", info.used)
quit()


# 500k jacobi
adata = sc.read_csv("/media/volume/sdb/mat500k.csv", delimiter=',', first_column_names= True, dtype='float32')

GPUtil.showUtilization()
from pynvml import *
nvmlInit()
handle = nvmlDeviceGetHandleByIndex(0)

start_time = time.time()
adata.obsm["X_pca"] = PCA(n_components=n_components, svd_solver = "jacobi", output_type="numpy").fit_transform(adata.X)
p = PCA(n_components=n_components, output_type="numpy")

sys.argv = time.time() - start_time
print("--- %s seconds ---" % int(sys.argv))
GPUtil.showUtilization()
info = nvmlDeviceGetMemoryInfo(handle)
print("Total memory:", info.total)
print("Free memory:", info.free)
print("Used memory:", info.used)
quit()



# 1M jacobi
adata = sc.read_csv("/media/volume/sdb/mat1M.csv", delimiter=',', first_column_names= True, dtype='float32')

GPUtil.showUtilization()
from pynvml import *
nvmlInit()
handle = nvmlDeviceGetHandleByIndex(0)

start_time = time.time()
adata.obsm["X_pca"] = PCA(n_components=n_components, svd_solver = "jacobi", output_type="numpy").fit_transform(adata.X)
p = PCA(n_components=n_components, output_type="numpy")

sys.argv = time.time() - start_time
print("--- %s seconds ---" % int(sys.argv))
GPUtil.showUtilization()
info = nvmlDeviceGetMemoryInfo(handle)
print("Total memory:", info.total)
print("Free memory:", info.free)
print("Used memory:", info.used)
quit()

# 1.3M jacobi
adata = sc.read_csv("/media/volume/sdb/mat13M.csv", delimiter=',', first_column_names= True, dtype='float32')

GPUtil.showUtilization()
from pynvml import *
nvmlInit()
handle = nvmlDeviceGetHandleByIndex(0)

start_time = time.time()
adata.obsm["X_pca"] = PCA(n_components=n_components, svd_solver = "jacobi", output_type="numpy").fit_transform(adata.X)
p = PCA(n_components=n_components, output_type="numpy")

sys.argv = time.time() - start_time
print("--- %s seconds ---" % int(sys.argv))
GPUtil.showUtilization()
info = nvmlDeviceGetMemoryInfo(handle)
print("Total memory:", info.total)
print("Free memory:", info.free)
print("Used memory:", info.used)
quit()



