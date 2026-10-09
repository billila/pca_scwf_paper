# How to extend the PCA benchmark

Guidelines for adding a **new dataset** or a **new PCA/SVD method** to the benchmark in [`pca/`](../pca).

---

## Adding a new dataset

1. **Choose the dataset**
   - Any count matrix (genes × cells) that can be loaded as a `SingleCellExperiment` (R) or `AnnData` (Python).
   - Prefer large datasets: the benchmark targets scalability.

2. **Preprocessing**
   - Log-normalize the counts (e.g. `scuttle::logNormCounts()`), as in [`pca/preprocessing/preprocessing_R_tenx.R`](../pca/preprocessing/preprocessing_R_tenx.R).
   - Choose the feature set:
     - select highly variable genes (HVG), **or**
     - keep all genes (this tests the methods on a wider matrix).
   - If the matrix is already normalized / log-transformed, skip this step and document it (see [`preprocessing_4.2M_cells.md`](../pca/preprocessing/preprocessing_4.2M_cells.md)).

3. **Subsampling for scalability (optional)**
   - Decide whether to create subsets of increasing size (e.g. 100k, 500k, 1M cells, full dataset).
   - Fix the random seed (`set.seed()`) so the same cells are used by every method.

4. **Save the input matrices**
   - Follow the scripts in [`pca/preprocessing/`](../pca/preprocessing) to save each subset in all the representations tested:
     - in-memory **dense** → `save_dense_matrix.R`
     - in-memory **sparse** → `save_sparse_matrix.R`
     - on-disk **HDF5 dense** → `save_hdf5_dense.R`
     - on-disk **HDF5 sparse** → `save_hdf5_sparse.R`
   - For the Python methods reading CSV (scikit-learn IPCA, RAPIDS), also export the subsets as `.csv` (cells × genes, float32).

5. **Run the PCA for all methods**
   - Update the input paths in the scripts of [`pca/run_pca_time/`](../pca/run_pca_time) and [`pca/run_pca_mem/`](../pca/run_pca_mem).
   - Run every R and Python implementation, **one at a time**, on a dedicated node.
   - Time: measured inside each script, around the PCA call only.
   - Memory: wrap each run with `/usr/bin/time -v` and record the `Maximum resident set size`.

6. **Collect the results**
   - Store elapsed time and peak memory per method × subset in a table with the same format as the one used for the paper figures ([`paper_figure/`](../paper_figure)).

---

## Adding a new PCA method

1. **Create one script per matrix representation** you want to test, named `<library>_<representation>_<algorithm>.{R,py}` (e.g. `newlib_sparse_randomized.py`).

2. **Use the same settings as the existing methods**
   - same input subsets (see above);
   - **50** principal components;
   - centering on, no scaling.

3. **Measure time and memory in the same way**
   - time only the PCA call (exclude data loading);
   - add a copy of the script to `run_pca_mem/` and run it with `/usr/bin/time -v`.

4. **Check the accuracy**
   - compare the PCs with those of an exact implementation (e.g. `BiocSingular::ExactParam()` or `sklearn PCA(svd_solver="full")`) on the same subset.

5. **Document the method**
   - add it to the implementation table in [`pca/README.md`](../pca/README.md);
   - add any new dependency to the environment files in [`envs/`](../envs).