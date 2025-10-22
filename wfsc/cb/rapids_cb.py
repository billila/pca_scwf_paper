# library 
import numpy as np
import pandas as pd 
import scanpy as sc
import time
import anndata
import os
import cudf
import cupy as cp
from sklearn.metrics import adjusted_rand_score
from sklearn.metrics import silhouette_score
from cuml.decomposition import PCA
from cuml.manifold import TSNE
from cuml.cluster import KMeans
from cuml.preprocessing import StandardScaler
import GPUtil
import rapids_scanpy_funcs
from scipy import sparse



sc.settings.verbosity = 3             # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_header()
sc.settings.set_figure_params(dpi=80, facecolor='white')

# save time usage #### 
time_sc = pd.DataFrame(index=["find_mit_gene", "filter", "normalization", "hvg",
                              "scaling", "PCA", "t-sne", "umap", "louvain", "leiden"],
                       columns=["time_sec"])



# data ####

adata = sc.read_h5ad("cord_blood.h5ad")
adata.var_names_make_unique()  
adata

genes = cudf.Series(adata.var_names)
adata.X = adata.X.astype('float32')
X_sparse = sparse.csr_matrix(adata.X)
celltype = cudf.Series(adata.obs['celltype'])
cells = cudf.Series(adata.obs_names)
X_sparse.colnames = adata.obs_names
print("X_sparse.colnames", X_sparse.colnames)
sparse_gpu_array = cp.sparse.csr_matrix(X_sparse)
sparse_gpu_array.colnames = adata.obs_names
print("sparse_gpu_array.colnames", sparse_gpu_array.colnames)
sparse_gpu_array.shape
print(sparse_gpu_array.shape)

# find mitocondrial genes ####
start_time = time.time()



end_time = time.time()
time_elapsed = end_time - start_time
print("Time Elapsed:", time_elapsed)
time_sc.iloc[0, 0] = time_elapsed

# filter data ####
start_time = time.time()

sparse_gpu_array, cells = rapids_scanpy_funcs.filter_cells_2(sparse_gpu_array, cells_idx=cells, min_genes=200, max_genes=2500)
print("cell dim:",cells.shape)
print("count matrix dim:", sparse_gpu_array.shape)

sparse_gpu_array, genes = rapids_scanpy_funcs.filter_genes(sparse_gpu_array, genes, min_cells=10)
print("count matrix dm 2:", sparse_gpu_array.shape)

end_time = time.time()
time_elapsed = end_time - start_time
print("Time Elapsed:", time_elapsed)
time_sc.iloc[1, 0] = time_elapsed

# normalization ####
start_time = time.time()

sparse_gpu_array = rapids_scanpy_funcs.normalize_total(sparse_gpu_array, target_sum=1e4)
# log transform
sparse_gpu_array = sparse_gpu_array.log1p()

end_time = time.time()
time_elapsed = end_time - start_time
print("Time Elapsed:", time_elapsed)
time_sc.iloc[2, 0] = time_elapsed

# Identification of highly variable features (feature selection) ####
start_time = time.time()



hvg = rapids_scanpy_funcs.highly_variable_genes(sparse_gpu_array, genes, n_top_genes=1000)


sparse_gpu_array = sparse_gpu_array[:, hvg]
genes = genes[hvg].reset_index(drop=True)
sparse_gpu_array.shape
sparse_gpu_array
end_time = time.time()
time_elapsed = end_time - start_time
print("Time Elapsed:", time_elapsed)
time_sc.iloc[3, 0] = time_elapsed
print("sparse_gpu_array.colnames", sparse_gpu_array.colnames)

# Scaling the data ####
start_time = time.time()

# sparse_gpu_array = cp.clip(StandardScaler().fit_transform(sparse_gpu_array), a_min = -10, a_max=10, with_mean=False)
print("sparse_gpu_array.colnames", sparse_gpu_array.colnames)

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


print("sparse_gpu_array.colnames", sparse_gpu_array.colnames)

adata = anndata.AnnData(sparse_gpu_array.get())
adata.var_names = genes.to_pandas()
adata.obs_names = cells.to_pandas()

#common_elements = adata.obs_names.intersection(celltype.index)
common_elements = sparse_gpu_array.obs_names.intersection(celltype.index)

print("dim common elements:", common_elements.shape)
new_celltype = celltype.loc[common_elements]
print("dim new celltype:", new_celltype.shape)

adata.obs['celltype'] = new_celltype



# PCA ####
start_time = time.time()
adata.obsm["X_pca"] = PCA(n_components=50, output_type="numpy").fit_transform(adata.X)
end_time = time.time()
time_elapsed = end_time - start_time
print("Time Elapsed:", time_elapsed)
time_sc.iloc[5, 0] = time_elapsed

# t-sne ####
start_time = time.time()
tsne_n_pcs = 20
adata.obsm['X_tsne'] = TSNE().fit_transform(adata.obsm["X_pca"][:,:tsne_n_pcs])

end_time = time.time()
time_elapsed = end_time - start_time
print("Time Elapsed:", time_elapsed)
time_sc.iloc[6, 0] = time_elapsed

# UMAP ####
start_time = time.time()
n_neighbors = 15 # Number of nearest neighbors for KNN graph
knn_n_pcs = 50 # Number of principal components to use for finding nearest neighbors
sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=knn_n_pcs, method='rapids')
umap_min_dist=0.3
umap_spread=1.0
sc.tl.umap(adata, min_dist=umap_min_dist, spread=umap_spread, method='rapids')

end_time = time.time()
time_elapsed = end_time - start_time
print("Time Elapsed:", time_elapsed)
time_sc.iloc[7, 0] = time_elapsed

# louvain ####
start_time = time.time()

n_neighbors=25
knn_n_pcs=50
sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=knn_n_pcs, method='rapids')
sc.tl.louvain(adata, flavor='rapids', resolution = 0.15)

end_time = time.time()
time_elapsed = end_time - start_time
print("Time Elapsed:", time_elapsed)
time_sc.iloc[8, 0] = time_elapsed

true_labels = adata.obs['celltype'].astype(str)
predicted_labels = adata.obs['louvain'].astype(str)

# Compute the ARI
ari_score = adjusted_rand_score(true_labels, predicted_labels)

# Print the ARI score
print("Adjusted Rand Index (ARI):", ari_score)

cluster_labels = adata.obs['louvain']

# Calculate silhouette scores
silhouette_avg = silhouette_score(adata.X, cluster_labels)

# Print the average silhouette score
print("Average Silhouette Score:", silhouette_avg)


# leiden ####
start_time = time.time()

adata.obs['leiden'] = rapids_scanpy_funcs.leiden(adata, resolution = 0.15)

end_time = time.time()
time_elapsed = end_time - start_time
print("Time Elapsed:", time_elapsed)
time_sc.iloc[9, 0] = time_elapsed

true_labels = adata.obs['celltype'].astype(str)
predicted_labels = adata.obs['leiden'].astype(str)

# Compute the ARI
ari_score = adjusted_rand_score(true_labels, predicted_labels)

# Print the ARI score
print("Adjusted Rand Index (ARI):", ari_score)

cluster_labels = adata.obs['leiden']

# Calculate silhouette scores
silhouette_avg = silhouette_score(adata.X, cluster_labels)

# Print the average silhouette score
print("Average Silhouette Score:", silhouette_avg)

print(time_sc)
time_sc
