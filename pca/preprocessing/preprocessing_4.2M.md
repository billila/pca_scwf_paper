# Preprocessing — 4.2M-cell MERFISH dataset
 
For the 4.2M-cell dataset (Zhuang-ABCA-1, 1,122 genes, [Allen Brain Cell Atlas](https://alleninstitute.github.io/abc_atlas_access/)) we used the **expression matrix already provided in log scale** by the Allen Institute:
 
```
Zhuang-ABCA-1-log2.h5ad
```
 
downloaded from:
https://allen-brain-cell-atlas.s3.us-west-2.amazonaws.com/index.html#expression_matrices/Zhuang-ABCA-1/20230830/
 
**No additional preprocessing was performed** on this matrix (no further filtering, normalization or log transformation): the log-transformed matrix was used directly as input for the analyses.
 
See also [`wfsc/input_data/4.2M_input.R`](../../wfsc/input_data/4.2M_input.R).