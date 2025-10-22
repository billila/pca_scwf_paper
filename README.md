# Benchmarking large-scale single-cell RNA-seq analysis

This repository contains the scripts and resources used in our paper  
**"Benchmarking large-scale single-cell RNA-seq analysis"**.  
It provides all materials necessary to reproduce the benchmarking results presented in the manuscript.

---

##  Reproducibility and environment setup

To ensure full reproducibility, we recommend recreating the computational environments used in our analyses.

- **R environment**  
  All R analyses were run within a pre-built container, available [here](https://github.com/billila/pca_scwf_paper/blob/main/envs/bioc_3_20_pca_wfsc.def).

- **Python environment**  
  For Python-based methods, we used a conda environment.  
  The corresponding `spca.yml` file can be found [here](https://github.com/billila/pca_scwf_paper/blob/main/envs/spca.yml).
  To create it:  

```
conda env create -f spca.yml
conda activate spca
```

---

## Repository Structure

```
├── pca/
│ ├── preprocessing/ # Scripts for data preprocessing prior to PCA
│ ├── run_pca/ # Code for all 28 PCA implementations benchmarked
│ └── README.md # Additional details on PCA benchmarking
│
├── scwf/
│ ├── 1.3M/ # Workflow scripts grouped by input dataset
│ ├── BE1/
│ ├── cb/
│ ├── sc_mix/
│ └── README.md
│
└── envs/ # Container and conda environment definitions
└── paper_figure/ #Code to reproduce paper figure
```


## Dataset

The datasets used in this benchmark are publicly available:

- **1.3M Brain Cells** — available as part of the [TENxBrainData](https://bioconductor.org/packages/TENxBrainData) Bioconductor package or 
you can register and download it from here: [1M_neurons.h5](https://support.10xgenomics.com/single-cell-gene-expression/datasets/1.3.0/1M_neurons)  
- **sc_mixology dataset** — available at [GEO: GSE118767](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE118767);  
  additional technical details are provided at [https://github.com/LuyiTian/sc_mixology](https://github.com/LuyiTian/sc_mixology).  
- **BE1 dataset** — available at [Figshare DOI: 10.6084/m9.figshare.23939481.v1](https://doi.org/10.6084/m9.figshare.23939481.v1).  
- **Cord Blood dataset** — available as part of the [SingleCellMultiModal](https://bioconductor.org/packages/SingleCellMultiModal) Bioconductor package.

---

## Citation

If you use this repository, please cite our paper:

> Billato *et al.* (2025). *Benchmarking large-scale single-cell RNA-seq analysis.*  
> [Journal / preprint link to be added]

---


## Session Info

<details>
<summary>
Click here for Session Info
</summary>

``` r
sessionInfo()
# R version 4.4.2 (2024-10-31)
# Platform: x86_64-pc-linux-gnu
# Running under: Ubuntu 24.04.1 LTS
# 
# Matrix products: default
# BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
# LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
# 
# locale:
#  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
#  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
#  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
#  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
#  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
# [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
# 
# time zone: Etc/UTC
# tzcode source: system (glibc)
# 
# attached base packages:
# [1] stats4    stats     graphics  grDevices utils     datasets  methods  
# [8] base     
# 
# other attached packages:
#  [1] tidyr_1.3.1                 patchwork_1.3.0            
#  [3] Seurat_5.2.1                SeuratObject_5.0.2         
#  [5] sp_2.2-0                    dplyr_1.1.4                
#  [7] scrapper_1.0.3              bluster_1.16.0             
#  [9] mclust_6.1.1                AnnotationDbi_1.68.0       
# [11] SingleCellMultiModal_1.18.0 MultiAssayExperiment_1.32.0
# [13] TENxBrainData_1.26.0        RSpectra_0.16-2            
# [15] rARPACK_0.11-0              DelayedMatrixStats_1.28.1  
# [17] BiocParallel_1.40.0         scran_1.34.0               
# [19] scater_1.34.0               scuttle_1.16.0             
# [21] ggplot2_3.5.1               mbkmeans_1.22.0            
# [23] HDF5Array_1.34.0            rhdf5_2.50.2               
# [25] DelayedArray_0.32.0         SparseArray_1.6.2          
# [27] S4Arrays_1.6.0              abind_1.4-8                
# [29] Matrix_1.7-2                here_1.0.1                 
# [31] BiocSingular_1.22.0         zellkonverter_1.16.0       
# [33] SingleCellExperiment_1.28.1 SummarizedExperiment_1.36.0
# [35] Biobase_2.66.0              GenomicRanges_1.58.0       
# [37] GenomeInfoDb_1.42.3         IRanges_2.40.1             
# [39] S4Vectors_0.44.0            BiocGenerics_0.52.0        
# [41] MatrixGenerics_1.18.1       matrixStats_1.5.0          
# 
# loaded via a namespace (and not attached):
#   [1] spatstat.sparse_3.1-0    httr_1.4.7               RColorBrewer_1.1-3      
#   [4] doParallel_1.0.17        tools_4.4.2              sctransform_0.4.1       
#   [7] R6_2.6.1                 lazyeval_0.2.2           uwot_0.2.3              
#  [10] rhdf5filters_1.18.0      withr_3.0.2              gridExtra_2.3           
#  [13] progressr_0.15.1         cli_3.6.4                spatstat.explore_3.3-4  
#  [16] fastDummies_1.7.5        spatstat.data_3.1-4      ggridges_0.5.6          
#  [19] pbapply_1.7-2            parallelly_1.42.0        limma_3.62.2            
#  [22] RSQLite_2.3.9            generics_0.1.3           ica_1.0-3               
#  [25] spatstat.random_3.3-2    ggbeeswarm_0.7.2         lifecycle_1.0.4         
#  [28] yaml_2.3.10              edgeR_4.4.2              BiocFileCache_2.14.0    
#  [31] Rtsne_0.17               grid_4.4.2               blob_1.2.4              
#  [34] promises_1.3.2           dqrng_0.4.1              ExperimentHub_2.14.0    
#  [37] crayon_1.5.3             dir.expiry_1.14.0        miniUI_0.1.1.1          
#  [40] lattice_0.22-6           beachmat_2.22.0          cowplot_1.1.3           
#  [43] KEGGREST_1.46.0          magick_2.8.5             pillar_1.10.1           
#  [46] metapod_1.14.0           rjson_0.2.23             future.apply_1.11.3     
#  [49] codetools_0.2-20         glue_1.8.0               spatstat.univar_3.1-1   
#  [52] data.table_1.17.0        vctrs_0.6.5              png_0.1-8               
#  [55] spam_2.11-1              gtable_0.3.6             cachem_1.1.0            
#  [58] mime_0.12                survival_3.8-3           iterators_1.0.14        
#  [61] statmod_1.5.0            gmp_0.7-5                fitdistrplus_1.2-2      
#  [64] ROCR_1.0-11              nlme_3.1-167             bit64_4.6.0-1           
#  [67] filelock_1.0.3           RcppAnnoy_0.0.22         rprojroot_2.0.4         
#  [70] irlba_2.3.5.1            vipor_0.4.7              KernSmooth_2.23-26      
#  [73] colorspace_2.1-1         DBI_1.2.3                tidyselect_1.2.1        
#  [76] bit_4.5.0.1              compiler_4.4.2           curl_6.2.1              
#  [79] BiocNeighbors_2.0.1      basilisk.utils_1.18.0    plotly_4.10.4           
#  [82] scales_1.3.0             lmtest_0.9-40            rappdirs_0.3.3          
#  [85] stringr_1.5.1            SpatialExperiment_1.16.0 digest_0.6.37           
#  [88] goftest_1.2-3            spatstat.utils_3.1-2     benchmarkmeData_1.0.4   
#  [91] basilisk_1.18.0          XVector_0.46.0           htmltools_0.5.8.1       
#  [94] pkgconfig_2.0.3          sparseMatrixStats_1.18.0 dbplyr_2.5.0            
#  [97] fastmap_1.2.0            rlang_1.1.5              htmlwidgets_1.6.4       
# [100] UCSC.utils_1.2.0         shiny_1.10.0             farver_2.1.2            
# [103] zoo_1.8-13               jsonlite_1.9.0           magrittr_2.0.3          
# [106] GenomeInfoDbData_1.2.13  dotCall64_1.2            Rhdf5lib_1.28.0         
# [109] munsell_0.5.1            Rcpp_1.0.14              viridis_0.6.5           
# [112] reticulate_1.41.0        stringi_1.8.4            ClusterR_1.3.3          
# [115] zlibbioc_1.52.0          MASS_7.3-64              AnnotationHub_3.14.0    
# [118] plyr_1.8.9               parallel_4.4.2           listenv_0.9.1           
# [121] ggrepel_0.9.6            deldir_2.0-4             Biostrings_2.74.1       
# [124] splines_4.4.2            tensor_1.5               locfit_1.5-9.11         
# [127] igraph_2.1.4             spatstat.geom_3.3-5      RcppHNSW_0.6.0          
# [130] reshape2_1.4.4           ScaledMatrix_1.14.0      BiocVersion_3.20.0      
# [133] BiocManager_1.30.25      foreach_1.5.2            httpuv_1.6.15           
# [136] RANN_2.6.2               purrr_1.0.4              polyclip_1.10-7         
# [139] future_1.34.0            benchmarkme_1.0.8        scattermore_1.2         
# [142] rsvd_1.0.5               xtable_1.8-4             later_1.4.1             
# [145] viridisLite_0.4.2        tibble_3.2.1             memoise_2.0.1           
# [148] beeswarm_0.4.0           cluster_2.1.8            globals_0.16.3 
```

</details>



