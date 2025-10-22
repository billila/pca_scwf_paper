import numpy as np
import scanpy as sc
import anndata
import pandas as pd
import time
import os
import cupy as cp
from cuml.decomposition import PCA
from cuml.manifold import TSNE
from cuml.cluster import KMeans
from cuml.preprocessing import StandardScaler

import rapids_scanpy_funcs

import warnings
warnings.filterwarnings('ignore', 'Expected ')
warnings.simplefilter('ignore')

import scipy
print(scipy.__version__)

import rmm

rmm.reinitialize(managed_memory=True)
cp.cuda.set_allocator(rmm.rmm_cupy_allocator)

input_file = "1M_brain_cells_10X.sparse.h5ad"

# maximum number of cells to load from files
USE_FIRST_N_CELLS = 1300000

# marker genes
MITO_GENE_PREFIX = "mt-"


# filtering cells
min_genes_per_cell = 200 
max_genes_per_cell = 6000 

# filtering genes
min_cells_per_gene = 1 
n_top_genes = 1000 

# PCA
n_components = 50 

# Batched PCA
pca_train_ratio = 0.35 
n_pca_batches = 10

# t-SNE
tsne_n_pcs = 20 

# k-means
k = 35 

# KNN
n_neighbors = 15 
knn_n_pcs = 50 

# UMAP
umap_min_dist = 0.3 
umap_spread = 1.0

############################
time_sc = pd.DataFrame(index=["find_mit_gene", "filter", "normalization", "hvg",
                              "scaling", "PCA", "t-sne", "umap", "louvain", "leiden"],
                       columns=["time_sec"])

###############################

# load data and filter gene and cells ##########
start_time = time.time()
sparse_gpu_array, genes, marker_genes_raw = \
    rapids_scanpy_funcs.preprocess_in_batches(input_file, 
                                              markers, 
                                              min_genes_per_cell=min_genes_per_cell, 
                                              max_genes_per_cell=max_genes_per_cell, 
                                              min_cells_per_gene=min_cells_per_gene, 
                                              target_sum=1e4, 
                                              n_top_genes=n_top_genes,
                                              max_cells=USE_FIRST_N_CELLS)
end_time = time.time()
time_elapsed = end_time - start_time
print("Time Elapsed:", time_elapsed)
time_sc.iloc[1, 0] = time_elapsed

print("step1",sparse_gpu_array)

# filter mithocondrial genes  

start_time = time.time()

mito_genes = genes.str.startswith(MITO_GENE_PREFIX)
print("stepmito_gene",sparse_gpu_array)

n_counts = sparse_gpu_array.sum(axis=1)
print("step_n_counts",sparse_gpu_array)

percent_mito = (sparse_gpu_array[:,mito_genes].sum(axis=1) / n_counts).ravel()
print("step_percent_mit",sparse_gpu_array)

n_counts = cp.array(n_counts).ravel()
print("step_ncount2",sparse_gpu_array)

percent_mito = cp.array(percent_mito).ravel()
print("step_percent_mito2",sparse_gpu_array)

end_time = time.time()
time_elapsed = end_time - start_time
print("Time Elapsed:", time_elapsed)
time_sc.iloc[3, 0] = time_elapsed  

# regress out #########

#sparse_gpu_array = rapids_scanpy_funcs.regress_out(sparse_gpu_array.tocsc(), n_counts, percent_mito)
#del n_counts, percent_mito, mito_genes

print("regress out", sparse_gpu_array)
# scaling 
start_time = time.time()

mean = sparse_gpu_array.mean(axis=0)
sparse_gpu_array -= mean
stddev = cp.sqrt(sparse_gpu_array.var(axis=0))
sparse_gpu_array /= stddev
sparse_gpu_array = cp.clip(sparse_gpu_array, a_min = -10, a_max=10)
del mean, stddev

end_time = time.time()
time_elapsed = end_time - start_time
print("Time Elapsed:", time_elapsed)
time_sc.iloc[4, 0] = time_elapsed  
print("scaling", sparse_gpu_array)
# new anndata
adata = anndata.AnnData(sparse_gpu_array.get())
print("adata", adata.X)
adata.var_names = genes.to_pandas()
del sparse_gpu_array

for name, data in marker_genes_raw.items():
    adata.obs[name] = data.get()
    
print(adata)
print(adata.X)
# PCA
 
start_time = time.time()

adata.obsm["X_pca"] = PCA(n_components=n_components, output_type="numpy").fit_transform(adata.X)

end_time = time.time()
time_elapsed = end_time - start_time
print("Time Elapsed:", time_elapsed)
time_sc.iloc[5, 0] = time_elapsed 

# t-sne
start_time = time.time()

adata.obsm['X_tsne'] = TSNE().fit_transform(adata.obsm["X_pca"][:,:tsne_n_pcs])

end_time = time.time()
time_elapsed = end_time - start_time
print("Time Elapsed:", time_elapsed)
time_sc.iloc[6, 0] = time_elapsed 

# umap 
sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=knn_n_pcs, method='rapids')

sc.tl.umap(adata, min_dist=umap_min_dist, spread=umap_spread, method='rapids')

end_time = time.time()
time_elapsed = end_time - start_time
print("Time Elapsed:", time_elapsed)
time_sc.iloc[7, 0] = time_elapsed

# louvain 

sc.tl.louvain(adata, flavor='rapids')
end_time = time.time()
time_elapsed = end_time - start_time
print("Time Elapsed:", time_elapsed)
time_sc.iloc[8, 0] = time_elapsed

# leiden 

adata.obs['leiden'] = rapids_scanpy_funcs.leiden(adata)

end_time = time.time()
time_elapsed = end_time - start_time
print("Time Elapsed:", time_elapsed)
time_sc.iloc[9, 0] = time_elapsed


print(time_sc)
time_sc
