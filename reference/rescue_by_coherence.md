# Rescue borderline/removed cells by cluster coherence

Rescue borderline/removed cells by cluster coherence

## Usage

``` r
rescue_by_coherence(
  obj,
  metrics,
  detected,
  celltype_col = "auto_celltype",
  min_cluster_size = 50,
  rescue_mode = c("none", "lenient", "moderate", "strict"),
  cancer_bypass = FALSE,
  coherence_min = NULL,
  gates_min_pass = NULL
)
```

## Arguments

- obj:

  Seurat object

- metrics:

  Character vector of metric names to consider.

- detected:

  Numeric vector of detected gene counts per cell (length = ncol(obj)).

- celltype_col:

  optional meta column with labels (e.g. 'celltype'); if missing,
  label-coherence is skipped

- min_cluster_size:

  minimum cluster size to consider for rescue (default 25)

- rescue_mode:

  "moderate" (default), "conservative", or "aggressive"

- coherence_min:

  Minimum coherence score required to rescue (default 0.5).

- gates_min_pass:

  Minimum number of gates a cell must pass to be rescued (default 1).

## Value

same `base` with possibly updated `qc_status`; adds `rescue_flag` and
`rescue_reason`
