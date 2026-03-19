# FeaturePlot-like map of auto-annotation confidence

FeaturePlot-like map of auto-annotation confidence

## Usage

``` r
scq_plot_auto_confidence(
  obj,
  imgDir = "qc_outputs",
  reduction = NULL,
  col = NULL,
  filename = "auto_confidence.png"
)
```

## Arguments

- obj:

  Seurat object

- imgDir:

  directory to write the PNG (created if missing)

- reduction:

  which embedding to use; default: first of umap/tsne/pca found

- col:

  meta column name for confidence; default tries common names

- filename:

  output filename

## Value

full path to the PNG (character) or NULL if skipped
