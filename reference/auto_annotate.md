# Auto-annotate cell types (multi-backend wrapper)

Auto-annotate cell types (multi-backend wrapper)

## Usage

``` r
auto_annotate(
  obj,
  species = c("mouse", "human"),
  assay = "RNA",
  prefer = c("sctype", "marker_score", "singleR", "label_transfer"),
  level = c("cluster", "cell"),
  marker_score_level = c("cluster", "cell"),
  custom_signatures = NULL,
  ref_seurat = NULL,
  ref_label_col = NULL,
  unknown_min_genes = 2,
  unknown_min_margin = 0.01,
  unknown_min_top = 0.01,
  organ = NULL,
  sctype_db_path = NULL,
  marker_method = c("module_score", "findmarkers"),
  panglao_file = system.file("extdata", "PanglaoDB_markers_27_Mar_2020.tsv", package =
    utils::packageName()),
  panglao_tissue = NULL,
  panglao_ui_max = 0.2,
  panglao_canonical_only = TRUE,
  fm_top_n = 50,
  fm_padj_max = 0.05,
  overlap_mode = c("count", "jaccard")
)
```

## Arguments

- obj:

  Seurat object.

- species:

  Character, "mouse" or "human".

- assay:

  Assay to use (default "RNA").

- prefer:

  Optional character; preferred backend(s) to try first.

- level:

  "cell" or "cluster".

- marker_score_level:

  Level for marker_score ("cell"/"cluster").

- custom_signatures:

  Optional list/data.frame of custom markers.

- ref_seurat:

  Optional Seurat reference object (label transfer).

- ref_label_col:

  Column in `ref_seurat@meta.data` with labels.

- unknown_min_genes:

  Integer; minimum genes to avoid "unknown".

- unknown_min_margin:

  Numeric; min score margin to call a label.

- unknown_min_top:

  Integer; min top markers expressed to call a label.

- organ:

  Optional organ/tissue context string.

- sctype_db_path:

  Optional path to ScType flat DB (xlsx).

- marker_method:

  Marker backend when using marker-based annotation.

- panglao_file:

  Optional PanglaoDB CSV/XLSX file path.

- panglao_tissue:

  Optional tissue filter for PanglaoDB.

- panglao_ui_max:

  Integer; cap markers per type to avoid overplot.

- panglao_canonical_only:

  Logical; keep only canonical markers.

- fm_top_n:

  Integer; FindMarkers top N to keep per cluster.

- fm_padj_max:

  Numeric; padj cutoff for FindMarkers filtering.

- overlap_mode:

  How to combine/overlap multiple marker sources.

## Value

Seurat object with annotation columns added.
