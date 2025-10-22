import sklearn as sk
import h5py as h5
from sklearn.decomposition import IncrementalPCA
import time
import os
import sys




# 100k 
os.chdir("/mnt/spca/run_spca_2025/M1/time/metodo2/data/subset/TENxBrainDataSE/TENxBrainData_100k/")
dsname = "assay001"
fn = "assays.h5"
n_components = 50

start_time = time.time()

matref = h5.File(fn, mode = "r")

#chunk_size = 10
op = sk.decomposition.PCA(n_components = n_components, svd_solver = "full")


#batch_size = chunk_size
op.fit_transform(matref[dsname])
sys.argv = time.time() - start_time
print("--- %s seconds ---" % int(sys.argv))

# 500k 
os.chdir("/mnt/spca/run_spca_2025/M1/time/metodo2/data/subset/TENxBrainDataSE/TENxBrainData_500k/")
dsname = "assay001"
fn = "assays.h5"
n_components = 50

start_time = time.time()

matref = h5.File(fn, mode = "r")

#chunk_size = 10
op = sk.decomposition.PCA(n_components = n_components, svd_solver = "full")

#batch_size = chunk_size
op.fit_transform(matref[dsname])
sys.argv = time.time() - start_time
print("--- %s seconds ---" % int(sys.argv))

# 1M
os.chdir("/mnt/spca/run_spca_2025/M1/time/metodo2/data/subset/TENxBrainDataSE/TENxBrainData_1000k/")
dsname = "assay001"
fn = "assays.h5"
n_components = 50

start_time = time.time()

matref = h5.File(fn, mode = "r")

#chunk_size = 10
op = sk.decomposition.PCA(n_components = n_components, svd_solver = "full")

#batch_size = chunk_size
op.fit_transform(matref[dsname])
sys.argv = time.time() - start_time
print("--- %s seconds ---" % int(sys.argv))

# 1.3M
os.chdir("/mnt/spca/run_spca_2025/M1/time/metodo2/data/subset/TENxBrainDataSE/TENxBrainData_1.3M/")
dsname = "assay001"
fn = "assays.h5"
n_components = 50

start_time = time.time()

matref = h5.File(fn, mode = "r")

#chunk_size = 10
op = sk.decomposition.PCA(n_components = n_components, svd_solver = "full")

#batch_size = chunk_size
op.fit_transform(matref[dsname])
sys.argv = time.time() - start_time
print("--- %s seconds ---" % int(sys.argv))
