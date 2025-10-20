#' Subset a Seurat object by QC status (Seurat v5-safe)
#'
#' @param obj Seurat object
#' @param status_df data.frame or logical named vector
#'   - If data.frame: use 'qc_status' (keep/borderline/remove) or 'is_low_quality' (TRUE=bad)
#'   - If logical: TRUE = keep
#' @param keep_levels character vector of statuses to retain (default c("keep","borderline"))
#' @export
applyQCfilter <- function(obj, status_df, keep_levels = c("keep","borderline")) {
  stopifnot(inherits(obj, "Seurat"))

  if (is.logical(status_df) && is.null(dim(status_df))) {
    stopifnot(!is.null(names(status_df)))
    keep <- status_df
  } else if (is.data.frame(status_df)) {
    rn <- rownames(status_df)
    if (is.null(rn)) stop("status_df must have rownames = cell barcodes.")
    if ("qc_status" %in% colnames(status_df)) {
      keep <- status_df$qc_status %in% keep_levels
      names(keep) <- rn
    } else if ("is_low_quality" %in% colnames(status_df)) {
      keep <- !as.logical(status_df$is_low_quality)
      names(keep) <- rn
    } else {
      stop("applyQCfilter(): status_df must have 'qc_status' or 'is_low_quality' column.")
    }
  } else {
    stop("applyQCfilter(): status_df must be data.frame or logical named vector.")
  }

  # align to object cells
  keep <- keep[colnames(obj)]
  keep[is.na(keep)] <- FALSE

  if (!any(keep)) stop("applyQCfilter(): no cells to keep after filtering.")
  if (getOption("scQCenrich.debug", FALSE)) {
    message("[scQCenrich-debug] applyQCfilter(): keeping ", sum(keep), " of ", length(keep), " cells.")
  }

  obj <- subset(obj, cells = names(keep)[keep])
  return(obj)
}
