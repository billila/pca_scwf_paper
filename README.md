# Benchmarking Large-Scale Single-Cell RNA-seq Analysis

This repository contains the scripts and resources used in our paper  
**"Benchmarking large-scale single-cell RNA-seq analysis"**.  
It provides all materials necessary to reproduce the benchmarking results presented in the manuscript.

---

##  Reproducibility and environment setup

To ensure full reproducibility, we recommend recreating the computational environments used in our analyses.

- **R environment**  
  All R analyses were run within a pre-built container, available [here](LINK_TO_R_CONTAINER).

- **Python environment**  
  For Python-based methods, we used a conda environment.  
  The corresponding `environment.yml` file can be found [here](LINK_TO_CONDA_ENV).

---

## Repository Structure

```
├── pca/
│ ├── preprocessing/ # Scripts for data preprocessing prior to PCA
│ ├── methods/ # Code for all 28 PCA implementations benchmarked
│ └── README.md # Additional details on PCA benchmarking
│
├── scwf/
│ ├── <dataset_name>/ # Workflow scripts grouped by input dataset
│ └── README.md
│
└── envs/ # Container and conda environment definitions
```

## **Setup Instructions**  

### **Python Environment**  

For Python, use the Conda environment **`spca`**.  
To create it, use the `spca.yml` file:  

```
conda env create -f spca.yml
conda activate spca
```

### **R Environment**  
For R, use the def file `bioc_3_20_pca_wsc.def`  

## Dataset

The datasets used in this benchmark are publicly available:

- **1.3M Brain Cells** — available as part of the [TENxBrainData](https://bioconductor.org/packages/TENxBrainData) Bioconductor package.  
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





