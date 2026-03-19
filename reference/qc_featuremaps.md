# Make & save UMAP featuremaps for QC metrics

Make & save UMAP featuremaps for QC metrics

## Usage

``` r
qc_featuremaps(
  obj,
  metrics = c("nCount_QC", "nFeature_QC", "pctMT_QC", "MALAT1_frac_QC",
    "stress_score_QC", "intronic_frac_QC", "spliced_frac_QC", "unspliced_frac_QC",
    "u2s_ratio_QC", "dbl_score"),
  save_dir = "qc_outputs",
  prefix = "qc",
  ncol = 4,
  logify = c(nCount_QC = TRUE, nFeature_QC = TRUE, u2s_ratio_QC = TRUE, dbl_score =
    FALSE)
)
```

## Arguments

- obj:

  Seurat object (must have UMAP). Use .ensure_pca_umap() beforehand.

- metrics:

  Vector of meta column names to plot (only those present are used)

- save_dir:

  Where to write PNGs

- prefix:

  File prefix for images

- ncol:

  Number of columns in the grid image

- logify:

  Named logical vector: which metrics to log1p/ log10; defaults below
