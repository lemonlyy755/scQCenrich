# High-level wrapper that computes QC metrics, flags low-quality cells, (optionally) detects doublets and performs lightweight annotation, and renders an HTML report.

High-level wrapper that computes QC metrics, flags low-quality cells,
(optionally) detects doublets and performs lightweight annotation, and
renders an HTML report.

## Usage

``` r
run_qc_pipeline(
  obj,
  species = c("mouse", "human"),
  method = c("gmm", "threshold", "auto"),
  score_cutoff = 2,
  auto_k = 3,
  assay = "RNA",
  save_dir = "qc_outputs",
  report_file = "qc_report.html",
  spliced = NULL,
  unspliced = NULL,
  spliced_obj = NULL,
  spliced_assay = NULL,
  spliced_layer = "counts",
  unspliced_obj = NULL,
  unspliced_assay = NULL,
  unspliced_layer = "counts",
  doublets = c("none", "auto", "doubletfinder", "scdblfinder"),
  remove_doublets = TRUE,
  sample_col = NULL,
  annot_method = c("marker_score", "singleR", "sctype", "none"),
  marker_source = c("panglao", "internal", "custom"),
  panglao_file = NULL,
  canonical_only = TRUE,
  ui_max = 0.2,
  tissue = NULL,
  min_genes = 5,
  report_html = TRUE,
  theme_bootswatch = "minty",
  debug = FALSE,
  sctype_db_path = system.file("extdata", "ScTypeDB_full.xlsx", package = "scQCenrich",
    mustWork = FALSE),
  marker_method = c("findmarkers", "module_score"),
  annot_args = NULL,
  unknown_min_genes = 2,
  unknown_min_margin = 0.01,
  unknown_min_top = 0.01,
  qc_strength = c("auto", "default", "lenient", "strict"),
  rescue_mode = c("moderate", "lenient", "strict", "none"),
  cancer_bypass = FALSE,
  enrichment_plots = TRUE
)
```

## Arguments

- obj:

  Seurat object

- species:

  'mouse' or 'human'

- method:

  'threshold', 'gmm' or 'auto'

- score_cutoff:

  Integer cutoff used by 'threshold' mode

- auto_k:

  Ignored; kept for back-compat

- assay:

  Assay name

- save_dir:

  Directory for intermediate outputs

- report_file:

  Output HTML report path

- spliced, unspliced:

  Back-compat layer strings on the same object

- spliced_obj, unspliced_obj:

  Optional external Seurat objects

- spliced_assay, unspliced_assay:

  Assay names on external objects

- spliced_layer, unspliced_layer:

  Layer/slot on external objects

- doublets:

  'auto','doubletfinder','scdblfinder','none'

- remove_doublets:

  Logical; if TRUE, mark doublets to remove/borderline

- sample_col:

  Optional meta column with sample/group

- annot_method:

  'SingleR','marker_score','none'

- marker_source:

  'panglao', 'internal', 'custom' (only used when
  annot_method='marker_score')

- panglao_file:

  Optional path to Panglao TSV (defaults to inst/extdata if present)

- canonical_only, ui_max, min_genes:

  Panglao signature filtering knobs

- tissue:

  Character vector of target organs/tissues, e.g. c("Liver", "Kidney").

- report_html:

  Logical; if TRUE, render the report

- theme_bootswatch:

  Bootswatch theme for report

- debug:

  Logical for verbose messages

- sctype_db_path:

  Path to ScType Excel DB. Defaults to the packaged file

- marker_method:

  Either "module_score" (default) or "findmarkers".

- annot_args:

  List of extra args forwarded to
  [`auto_annotate()`](https://lemonlyy755.github.io/scQCenrich/reference/auto_annotate.md).

- unknown_min_genes:

  Integer; guard for "unknown" assignment.

- unknown_min_margin:

  Numeric; guard for "unknown" assignment.

- unknown_min_top:

  Integer; guard for "unknown" assignment.

- qc_strength:

  One of `c("auto","default","lenient","strict")`; passed to
  [`flagLowQuality()`](https://lemonlyy755.github.io/scQCenrich/reference/flagLowQuality.md).

- rescue_mode:

  One of `c("moderate","lenient","strict","none")`; controls how
  aggressively borderline cells are rescued. Default `"moderate"`.

- cancer_bypass:

  Logical. If TRUE, clusters with healthy splicing profiles but high
  removal rates are exempt from the removal-fraction penalty (useful for
  cancer datasets). Default `FALSE`.

- enrichment_plots:

  logical. If TRUE, run GO/KEGG enrichment plots; if FALSE, skip.
  Default: TRUE.

## Value

A list with at least: `$obj` (QC-kept), `$metrics`, `$status_df`,
`$report_file`. Also returns `$obj_all` for convenience.
