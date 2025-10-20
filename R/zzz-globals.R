#' @keywords internal
#' @importFrom methods as
#' @importFrom stats aggregate reorder setNames
#' @importFrom utils head
#' @importFrom grDevices png
#' @importFrom graphics axis image par
#' @importFrom SeuratObject DefaultAssay "DefaultAssay<-"
#' @importFrom rlang .data
#' @importFrom utils write.csv tail
#' @importFrom stats sd na.omit
#' @importFrom graphics plot.new text mtext



utils::globalVariables(c(
  ".", "variable", "value", "qc_status", "low_quality", "Freq",
  "pctMT_QC", "nFeature_QC", "stress_score_QC", "MALAT1_frac_QC", "intronic_frac_QC",
  "spliced_frac_QC", "unspliced_frac_QC", "u2s_ratio_QC",
  "gs_name", "gene_symbol", "cluster", "label", "UMAP_1","UMAP_2","val","is_doublet","group","Description","Count", "celltype","p_neglog10","term_label","gene",
  "x","y","conf","placeholder"
))


#' toy_seu
#'
#' A small example Seurat object used in examples and tests.
#'
#' @format A \code{Seurat} object with minimal counts and metadata.
#' @source Generated in package development for demos.
#' @usage data(toy_seu)
"toy_seu"

