# Removed cells clustering + enrichment (used by report)

Removed cells clustering + enrichment (used by report)

## Usage

``` r
removed_cells_analysis(
  obj,
  status_df,
  save_dir = "qc_outputs",
  species = c("mouse", "human"),
  assay = "RNA"
)
```

## Arguments

- obj:

  Seurat object

- status_df:

  Data frame with 'qc_status'

- save_dir:

  Output directory for PNG/CSV

- species:

  'mouse' or 'human'

- assay:

  Assay name

## Value

Invisibly returns NULL; writes files to disk
