# How to extend the single-cell workflow benchmark

Guidelines for adding a **new dataset** or a **new workflow** to the benchmark in [`wfsc/`](../wfsc).

---

## Adding a new dataset

1. **Identify the dataset**
   - Raw counts (genes × cells) are required, since preprocessing is performed inside each workflow.
   - Prefer datasets with **ground-truth cell labels** (cell lines, sorted populations, CITE-seq-based annotations), so clustering accuracy (ARI) can be computed.

2. **Create the input objects**
   - Add a script `wfsc/input_data/<dataset>_input.R`, following the existing ones, that:
     - downloads / loads the data;
     - saves a `SingleCellExperiment` (`.RData` / `.rds`) for the R workflows;
     - saves an `.h5ad` (`zellkonverter::writeH5AD()`) for the Python workflows.

3. **Create the dataset folder**
   - `wfsc/<dataset>/`, with one script per workflow: `OSCA_<dataset>.R`, `scrapper_<dataset>.R`, `seurat_<dataset>.R`, `scanpy_<dataset>.py`, `rapids_<dataset>.py`.
   - Start from the scripts of an existing dataset and change only the input file and the dataset-specific parameters.

4. **Preprocessing inside each workflow**
   - Each workflow performs its own QC, filtering, normalization, HVG selection (top **1,000** genes) and, where applicable, scaling, using the functions of its own package.
   - Adapt QC thresholds to the dataset (e.g. mitochondrial gene prefix `MT-` for human vs `mt-` for mouse).

5. **Run and collect the results**
   - Run one workflow at a time, wrapped with `/usr/bin/time -v` for peak memory.
   - Keep the per-step `time` table (10 steps: `find_mit_gene`, `filter`, `normalization`, `hvg`, `scaling`, `PCA`, `t-sne`, `umap`, `louvain`, `leiden`).
   - Compute the ARI between Louvain/Leiden clusters and the reference labels.

---

## Adding a new workflow

1. **Create one script per dataset**, named `<workflow>_<dataset>.{R,py}`, in each dataset folder.

2. **Implement the same steps** as the existing workflows, using the package's recommended functions:
   - QC and filtering → normalization → 1,000 HVG → (scaling) → PCA with **50** components → t-SNE → UMAP → Louvain → Leiden.
   - If a step is not part of the workflow, leave it as `NA` in the `time` table.

3. **Measure time and memory in the same way**
   - time each step separately (`Sys.time()` in R, `time.time()` in Python);
   - record peak memory with `/usr/bin/time -v`.

4. **Evaluate clustering**
   - compute the ARI against the reference labels (`mclust::adjustedRandIndex()` / `sklearn.metrics.adjusted_rand_score()`).

5. **Document the workflow**
   - add it to the workflow table in [`wfsc/README.md`](../wfsc/README.md);
   - add a container / conda environment file to [`envs/`](../envs).