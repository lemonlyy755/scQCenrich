# Subset a Seurat object by QC status (Seurat v5-safe)

Subset a Seurat object by QC status (Seurat v5-safe)

## Usage

``` r
applyQCfilter(obj, status_df, keep_levels = c("keep", "borderline"))
```

## Arguments

- obj:

  Seurat object

- status_df:

  data.frame or logical named vector

  - If data.frame: use 'qc_status' (keep/borderline/remove) or
    'is_low_quality' (TRUE=bad)

  - If logical: TRUE = keep

- keep_levels:

  character vector of statuses to retain (default
  c("keep","borderline"))
