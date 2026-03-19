# Ensure PCA/UMAP are present for nice plots (Seurat v5)

Ensure PCA/UMAP are present for nice plots (Seurat v5)

## Usage

``` r
.ensure_visual_embeddings(
  seu,
  assay = "RNA",
  pc_dims = 1:20,
  debug = getOption("scQCenrich.debug", FALSE)
)
```
