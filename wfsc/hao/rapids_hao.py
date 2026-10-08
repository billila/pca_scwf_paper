# hao rapids_singlecell workflow

import scanpy as sc
import cupy as cp
import pandas as pd
from sklearn.metrics import adjusted_rand_score
import time
import rapids_singlecell as rsc
import warnings
warnings.filterwarnings("ignore")

import rmm
from rmm.allocators.cupy import rmm_cupy_allocator
rmm.reinitialize(managed_memory=True, pool_allocator=False, devices=0)
cp.cuda.set_allocator(rmm_cupy_allocator)

time_sc = pd.DataFrame(index=["find_mit_gene", "filter", "normalization", "hvg",
                           "scaling", "PCA", "t-sne", "knn", "umap", "louvain", "leiden"],
                    columns=["time_sec"])

adata = sc.read_h5ad("cite_osca_filtered.h5ad")
print("Loaded CITE-seq:", adata.shape)
rsc.get.anndata_to_GPU(adata)
time_sc.iloc[0, 0] = 0; time_sc.iloc[1, 0] = 0

start_time = time.time()
rsc.pp.normalize_total(adata, target_sum=1e4)
rsc.pp.log1p(adata)
time_sc.iloc[2, 0] = time.time() - start_time

start_time = time.time()
adata.layers["counts"] = adata.X.copy()
rsc.pp.highly_variable_genes(adata, n_top_genes=1000, flavor="seurat_v3", layer="counts")
adata.raw = adata
adata = adata[:, adata.var["highly_variable"]]
time_sc.iloc[3, 0] = time.time() - start_time

start_time = time.time()
rsc.pp.scale(adata, max_value=10)
time_sc.iloc[4, 0] = time.time() - start_time

start_time = time.time()
rsc.pp.pca(adata, n_comps=50)
rsc.get.anndata_to_CPU(adata, convert_all=True)
time_sc.iloc[5, 0] = time.time() - start_time

start_time = time.time()
rsc.tl.tsne(adata, n_pcs=50)
time_sc.iloc[6, 0] = time.time() - start_time

start_time = time.time()
rsc.pp.neighbors(adata, n_neighbors=10, n_pcs=50)
time_sc.iloc[7, 0] = time.time() - start_time

start_time = time.time()
rsc.tl.umap(adata)
time_sc.iloc[8, 0] = time.time() - start_time

start_time = time.time()
rsc.tl.louvain(adata, resolution=0.5)
time_sc.iloc[9, 0] = time.time() - start_time
ari = adjusted_rand_score(adata.obs['celltype.l1'].astype(str), adata.obs['louvain'].astype(str))
print("Louvain ARI (celltype.l1):", ari)

start_time = time.time()
rsc.tl.leiden(adata, resolution=0.5)
time_sc.iloc[10, 0] = time.time() - start_time
ari = adjusted_rand_score(adata.obs['celltype.l1'].astype(str), adata.obs['leiden'].astype(str))
print("Leiden ARI (celltype.l1):", ari)

print(time_sc)
