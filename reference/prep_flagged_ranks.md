# Prepare ranked gene stats for dropped vs kept (logCPM mean diff)

Prepare ranked gene stats for dropped vs kept (logCPM mean diff)

## Usage

``` r
prep_flagged_ranks(obj, status_df, assay = "RNA")
```

## Arguments

- obj:

  Seurat object

- status_df:

  Data frame with rownames=cells and 'qc_status'

- assay:

  Assay name to use

## Value

Named numeric vector of statistics (higher = up in dropped)
