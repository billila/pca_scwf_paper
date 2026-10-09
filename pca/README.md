This folder contains the code used to benchmark the PCA/SVD implementations
available in R (Bioconductor, RSpectra) and Python (scanpy, scikit-learn,
RAPIDS/cuML) on the 1.3 Million Brain Cells 
dataset (10x Genomics, TENxBrainData) and the 4.2 MERFISH dataset.

Each implementation is run on four random subsets of increasing
size — 100k, 500k, 1M and 1.3M cells — and, for every subset,
we record elapsed time and peak memory usage. 
All methods compute the first 50 principal components on log-normalized
expression of the highly variable genes.

pca/

├── preprocessing/   # Gene filtering, normalization, subsetting and saving of the input matrices

├── run_pca_time/    # One script per implementation: elapsed time

└── run_pca_mem/     # One script per implementation: memory profiling

## Implementations benchmarked

Script names follow the pattern **`<library>_<matrix representation>_<algorithm>`**:

- **matrix representation**: `dense` (in-memory dense), `sparse` (in-memory sparse), `hdf5_dense` / `hdf5_sparse` (on-disk, `HDF5Array`/`DelayedArray`)
- **`_def`** suffix: `BiocSingular` with `deferred = TRUE`, i.e. centering is applied implicitly during matrix multiplication instead of creating a centered dense copy of the matrix


## Measuring time (`run_pca_time/`)

Each script loads the four subsets in turn and times **only the PCA call** (data loading excluded):

- R: `proc.time()` before/after `runPCA()` / `svds()`; elapsed time is reported in minutes.
- Python: `time.time()` before/after the PCA call; elapsed time is reported in seconds.
- RAPIDS: GPU memory usage is additionally printed via `pynvml`.

```bash
Rscript run_pca_time/bioc_sparse_irlba.R
python  run_pca_time/scanpy_sparse_arpack.py
```

## Measuring memory (`run_pca_mem/`)

The scripts in this folder mirror those in `run_pca_time/` (R implementations) and add memory profiling around the PCA call with `Rprof(memory.profiling = TRUE)`; profiles are written to `output/`.

Peak resident memory (max RSS) for **both R and Python** scripts was obtained by wrapping each run with GNU `time`, as shown in `save_mem_python_R.sh`:

```bash
/usr/bin/time -v -- Rscript run_pca_mem/bioc_sparse_irlba.R
/usr/bin/time -v -- python  run_pca_time/scanpy_sparse_arpack.py
```

and reading the `Maximum resident set size` line of the output.

## Notes

- Input paths are hard-coded to the directories used on our machines (e.g. `/mnt/spca/...`, `here("metodo2/...")`); adjust them to your local layout before running.
- Scripts are meant to be run one implementation at a time, on a dedicated node, so that time and memory measurements are not affected by other processes.

