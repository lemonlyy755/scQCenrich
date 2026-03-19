# Render a pretty HTML QC report

Render a pretty HTML QC report

## Usage

``` r
render_qc_report(
  obj,
  metrics,
  status_df,
  save_dir = "qc_outputs",
  report_file = file.path(save_dir, "qc_report.html"),
  theme_bootswatch = "minty",
  self_contained = TRUE,
  sample_col = NULL,
  template = NULL,
  species = c("mouse", "human")
)
```

## Arguments

- obj:

  Seurat object (with qc_status in meta)

- metrics:

  data.frame from calcQCmetrics()

- status_df:

  data.frame from flagLowQuality()

- save_dir:

  output dir for assets/plots (and where the Rmd will be copied)

- report_file:

  output HTML path (default: file.path(save_dir, "qc_report.html"))

- theme_bootswatch:

  Bootswatch theme name

- self_contained:

  embed all assets into single HTML (TRUE)

- sample_col:

  optional meta column for per-sample summaries

- template:

  path to a .Rmd template; if NULL, copy bundled template to save_dir
  and use it

- species:

  Character scalar, one of "mouse" or "human"
