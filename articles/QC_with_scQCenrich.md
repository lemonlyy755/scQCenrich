# QC_with_scQCenrich

## Overview

This vignette runs
[`run_qc_pipeline()`](https://lemonlyy755.github.io/scQCenrich/reference/run_qc_pipeline.md)
on Seurat objects shipped with the package under `inst/extdata/`.

``` r

library(scQCenrich)
library(Seurat)
```

``` r

std_path <- system.file("extdata","seurat_std_small.rds", package="scQCenrich")
sp_path  <- system.file("extdata","seurat_sp_small.rds", package="scQCenrich")
stopifnot(nzchar(std_path), file.exists(std_path))
stopifnot(nzchar(sp_path),  file.exists(sp_path))
seurat_std_small <- readRDS(std_path)
seurat_sp_small  <- readRDS(sp_path)
```

## Run the pipeline (fast, minimal settings)

``` r

res <- run_qc_pipeline(
  obj       = seurat_std_small,
  species   = "mouse",
  # external splicing objects:
  spliced_obj     = seurat_sp_small,
  spliced_assay   = "spliced",
  spliced_layer   = "counts",
  unspliced_obj   = seurat_sp_small,
  unspliced_assay = "unspliced",
  unspliced_layer = "counts",
  report_html = TRUE,
  report_file = "qc_outputs/qc_report.html",
  tissue = "Skin")
```

### Inspect outputs

``` r

# Key metadata columns
head( res$obj@meta.data[, c("auto_celltype","qc_status","is_doublet")], 10)

# If report_html=TRUE, a report is written to:
normalizePath("qc_outputs/qc_report.html", winslash = "/")
```
