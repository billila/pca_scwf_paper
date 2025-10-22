import scanpy as sc
import cupy as cp
import pandas as pd
import time
import rapids_singlecell as rsc
import warnings
warnings.filterwarnings("ignore")
import rmm
from rmm.allocators.cupy import rmm_cupy_allocator

rmm.reinitialize(
    managed_memory=True,  # Allows oversubscription
    pool_allocator=False,  # default is False
    devices=0,  # GPU device IDs to register. By default registers only GPU 0.
)

cp.cuda.set_allocator(rmm_cupy_allocator)

adata = sc.read("BE1.h5ad")
print(adata) 

# save time usage #### 
time_sc = pd.DataFrame(index=["find_mit_gene", "filter", "normalization", "hvg",
                           "scaling", "PCA", "t-sne", "umap", "louvain", "leiden"],
                    columns=["time_sec"])

rsc.get.anndata_to_GPU(adata)
print(adata.shape)

# find mitocondrial genes ####
start_time = time.time()

rsc.pp.flag_gene_family(adata, gene_family_name="MT", gene_family_prefix="MT-")

rsc.pp.flag_gene_family(adata, gene_family_name="RIBO", gene_family_prefix="RPS")

rsc.pp.calculate_qc_metrics(adata, qc_vars=["MT", "RIBO"])

# plot MT and RIBO
#sc.pl.scatter(adata, x="total_counts", y="pct_counts_MT")
#sc.pl.scatter(adata, x="total_counts", y="n_genes_by_counts")

end_time = time.time()
time_elapsed = end_time - start_time
print("Time Elapsed:", time_elapsed)
time_sc.iloc[0, 0] = time_elapsed

# filter 
start_time = time.time()

#sc.pl.violin(adata, "n_genes_by_counts", jitter=0.4, groupby="PatientNumber")
#sc.pl.violin(adata, "total_counts", jitter=0.4, groupby="PatientNumber")
#sc.pl.violin(adata, "pct_counts_MT", jitter=0.4, groupby="PatientNumber")
adata = adata[adata.obs["n_genes_by_counts"] > 200]
adata = adata[adata.obs["n_genes_by_counts"] < 5000]
print(adata.shape)

adata = adata[adata.obs["pct_counts_MT"] < 5]
print(adata.shape)

rsc.pp.filter_genes(adata, min_count=3)

adata.layers["counts"] = adata.X.copy()
print(adata.shape)

end_time = time.time()
time_elapsed = end_time - start_time
print("Time Elapsed:", time_elapsed)
time_sc.iloc[1, 0] = time_elapsed

# Normalization

start_time = time.time()

rsc.pp.normalize_total(adata, target_sum=1e4)
rsc.pp.log1p(adata)

end_time = time.time()
time_elapsed = end_time - start_time
print("Time Elapsed:", time_elapsed)
time_sc.iloc[2, 0] = time_elapsed

# HGV
start_time = time.time()

rsc.pp.highly_variable_genes(
    adata,
    n_top_genes=1000,
    flavor="seurat_v3",
    #batch_key="PatientNumber",
    layer="counts",
)
adata.raw = adata
adata = adata[:, adata.var["highly_variable"]]
adata.shape

rsc.pp.regress_out(adata, keys=["total_counts", "pct_counts_MT"])

end_time = time.time()
time_elapsed = end_time - start_time
print("Time Elapsed:", time_elapsed)
time_sc.iloc[3, 0] = time_elapsed

# scaling 
start_time = time.time()

rsc.pp.scale(adata, max_value=10)

end_time = time.time()
time_elapsed = end_time - start_time
print("Time Elapsed:", time_elapsed)
time_sc.iloc[4, 0] = time_elapsed

# PCA
start_time = time.time()
rsc.pp.pca(adata, n_comps=50)
# sc.pl.pca_variance_ratio(adata, log=True, n_pcs=50)
rsc.get.anndata_to_CPU(adata, convert_all=True)

end_time = time.time()
time_elapsed = end_time - start_time
print("Time Elapsed:", time_elapsed)
time_sc.iloc[5, 0] = time_elapsed

# t-sne
start_time = time.time()

rsc.tl.tsne(adata, n_pcs=40)


end_time = time.time()
time_elapsed = end_time - start_time
print("Time Elapsed:", time_elapsed)
time_sc.iloc[6, 0] = time_elapsed

# UMAP
start_time = time.time()

rsc.pp.neighbors(adata, n_neighbors=10, n_pcs=50)
rsc.tl.umap(adata)

end_time = time.time()
time_elapsed = end_time - start_time
print("Time Elapsed:", time_elapsed)
time_sc.iloc[7, 0] = time_elapsed

# louvain 
start_time = time.time()

rsc.tl.louvain(adata, resolution=0.6)

end_time = time.time()
time_elapsed = end_time - start_time
print("Time Elapsed:", time_elapsed)
time_sc.iloc[8, 0] = time_elapsed

# leiden 
start_time = time.time()

rsc.tl.leiden(adata, resolution=0.6)

end_time = time.time()
time_elapsed = end_time - start_time
print("Time Elapsed:", time_elapsed)
time_sc.iloc[9, 0] = time_elapsed


time_sc
print(time_sc)
