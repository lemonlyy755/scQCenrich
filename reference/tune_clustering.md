# Grid-search to preserve small clusters (k.param & resolution)

Grid-search to preserve small clusters (k.param & resolution)

## Usage

``` r
tune_clustering(
  obj,
  dims_use = 1:30,
  k_vals = c(10, 12, 15, 20),
  res_vals = seq(0.2, 1.6, by = 0.2),
  algorithm = 4,
  prune = 0.1,
  small_frac_target = c(0.005, 0.05),
  seed = 1
)
```

## Arguments

- obj:

  Seurat object

- dims_use:

  Integer vector of PC indices

- k_vals:

  Integer vector of k.param candidates

- res_vals:

  Numeric vector of resolution candidates

- algorithm:

  Integer; clustering algorithm (e.g., Leiden=4)

- prune:

  Numeric; prune.SNN

- small_frac_target:

  Length-2 numeric; desired min cluster fraction range

- seed:

  Integer random seed

## Value

List with components \$stats (grid) and \$pick (chosen row)
