# Apply chosen k/res to an object (returns updated object)

Apply chosen k/res to an object (returns updated object)

## Usage

``` r
apply_tuned_clustering(
  obj,
  dims_use,
  k_best,
  res_best,
  prune = 0.1,
  algorithm = 4,
  seed = 1
)
```

## Arguments

- obj:

  Seurat object to update

- dims_use:

  Integer vector of PC indices used to build the graph

- k_best:

  Integer; chosen k.param

- res_best:

  Numeric; chosen resolution

- prune:

  Numeric; prune.SNN passed to FindNeighbors

- algorithm:

  Integer; clustering algorithm (e.g., 4 = Leiden)

- seed:

  Integer random seed

## Value

Updated Seurat object
