# Single-cell workflow benchmark

This folder contains the code used to benchmark **complete single-cell
RNA-seq analysis workflows** implemented in R and Python,
from raw counts to clustering, on the five datesets. 

For every workflow and dataset we time each step of a standard analysis and,
where ground-truth labels are available, assess clustering accuracy
with the **Adjusted Rand Index (ARI)**.

---

## Folder structure

```
wfsc/
├── input_data/   # Scripts to download / build the input objects for each dataset
├── 1.3M/         # 1.3 Million Brain Cells (10x Genomics)
├── BE1/          # BE1 - lung cancer cell lines + PBMCs
├── cb/           # Cord Blood CITE-seq (RNA modality)
├── hao/          # Hao et al. PBMC CITE-seq (RNA modality)
└── sc_mix/       # sc_mixology - 5 cell lines (10x)
```

Each dataset folder contains one script per workflow, named `<workflow>_<dataset>.{R,py}`.

QC thresholds, HVG method and clustering parameters follow
each package's tutorial/recommended defaults; the resolution
used for Louvain/Leiden is set in each script.

Times are measured with `Sys.time()` (R) and `time.time()`
(Python) around each step. For datasets with annotated cell identities 
(sc_mix, BE1, cb, hao), the ARI between Louvain/Leiden clusters
and the reference labels is computed with `mclust::adjustedRandIndex()` 
(R) or `sklearn.metrics.adjusted_rand_score()` (Python).

Peak memory was recorded by wrapping each run with GNU `time`:

```bash
/usr/bin/time -v -- Rscript sc_mix/OSCA_sc_mix.R
```

## Notes

- Input and output paths are hard-coded to the directories used on our machines (e.g. `/mnt/spca/pipeline_sc/...`); adjust them before running.
- Diagnostic plots (`VlnPlot`, `sc.pl.*`) are included as in the original tutorials but are not part of the timed results.