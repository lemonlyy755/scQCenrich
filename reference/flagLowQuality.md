# Flag low-quality cells

Flag low-quality cells

## Usage

``` r
flagLowQuality(
  metrics,
  method = c("threshold", "gmm", "hdbscan"),
  score_cutoff = 2,
  auto_k = 3,
  obj = NULL,
  celltype_col = NULL,
  sample_col = NULL,
  rescue_mode = c("moderate", "lenient", "strict", "none"),
  min_cluster_size = NULL,
  doublet_action = c("remove", "borderline", "none"),
  doublet_col = "is_doublet",
  qc_strength = c("auto", "default", "lenient", "strict"),
  save_dir = NULL
)
```

## Arguments

- metrics:

  data.frame of QC metrics (rows=cells). If it contains a logical column
  'is_doublet' (as added by calcQCmetrics), it will be honored.

- method:

  "threshold", "gmm", or "hdbscan"

- score_cutoff:

  integer for threshold mode (default 2)

- auto_k:

  kept for API compat (unused here)

- obj:

  optional Seurat object for rescue step

- celltype_col:

  optional meta column for rescue label coherence

- sample_col:

  optional sample column (passed to detector; may be NULL)

- rescue_mode:

  "moderate","lenient","strict", or "none"

- min_cluster_size:

  minimum cluster size to consider for rescue; NULL = adaptive (~0.5% or
  \>=10)

- doublet_action:

  what to do with predicted doublets in `metrics$is_doublet`: "remove"
  (default), "borderline", or "none"

- doublet_col:

  name of the doublet flag column in `metrics` (default "is_doublet")

- qc_strength:

  One of c("auto","default","lenient","strict")

- save_dir:

  Character path or NULL. If non-NULL, write intermediate CSVs here.

## Value

data.frame with rownames=cells and columns: qc_status, qc_score,
cluster, reason,
