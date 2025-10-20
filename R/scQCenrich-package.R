#' scQCenrich: State-aware QC for single-cell RNA-seq (with spliced/unspliced support)
#'
#' Features:
#' * QC outlier detection (GMM/HDBSCAN/threshold) over mito%, genes, stress, MALAT1,
#'   and (optionally) spliced/unspliced metrics (intronic_frac, unspliced/spliced fractions & ratio).
#' * Rescue of coherent, non-extreme flagged cells (3-tier keep/borderline/drop).
#' * Doublet detection (DoubletFinder / scDblFinder).
#' * Automatic cell-type annotation
#' * Optional clustering auto-tuning to preserve small clusters (k.param & resolution).
#' * Presentable HTML report (QC×label, removed-cell UMAP, enrichment).
#'
#' See README for how to set `spliced` / `unspliced` inputs (assay or assay:layer).
#' @keywords internal
"_PACKAGE"
