# ---- R/rescue.R ----
# Safe rescue of good, coherent clusters; never writes zero-length columns.

#' Rescue borderline/removed cells by cluster coherence
#' @param obj Seurat object
#' @param celltype_col optional meta column with labels (e.g. 'celltype'); if missing, label-coherence is skipped
#' @param min_cluster_size minimum cluster size to consider for rescue (default 25)
#' @param rescue_mode "moderate" (default), "conservative", or "aggressive"
#' @param metrics Character vector of metric names to consider.
#' @param detected Numeric vector of detected gene counts per cell (length = ncol(obj)).
#' @param coherence_min Minimum coherence score required to rescue (default 0.5).
#' @param gates_min_pass Minimum number of gates a cell must pass to be rescued (default 1).
#' @return same `base` with possibly updated `qc_status`; adds `rescue_flag` and `rescue_reason`
#' @export
rescue_by_coherence <- function(
    obj, metrics, detected,
    celltype_col = "auto_celltype",
    min_cluster_size = 50,
    rescue_mode = c("none","lenient","moderate","strict"),
    # explicit gates (kept for compatibility, but now only used for debug table)
    coherence_min = NULL,       # e.g., 0.60 for moderate
    gates_min_pass = NULL       # how many metric gates (mito/feat/stress) must pass
){
  rescue_mode <- match.arg(rescue_mode)

  # ---------------- rescue score thresholds by mode -------------------------
  resc_threshold <- switch(rescue_mode,
                           "strict"   = 0.90,
                           "moderate" = 0.70,
                           "lenient"  = 0.50,
                           "none"     = Inf
  )

  # compatibility: derive legacy defaults for printing in the gate table
  mode_defaults <- switch(rescue_mode,
                          "strict"   = list(coh = 0.70, nmin = 3L),
                          "moderate" = list(coh = 0.60, nmin = 2L),
                          "lenient"  = list(coh = 0.50, nmin = 1L),
                          "none"     = list(coh = 1.00, nmin = 99L)
  )
  if (is.null(coherence_min))  coherence_min  <- mode_defaults$coh
  if (is.null(gates_min_pass)) gates_min_pass <- mode_defaults$nmin

  stopifnot(is.data.frame(detected), nrow(detected) > 0L)
  stopifnot(all(c("qc_status") %in% colnames(detected)))

  # ---- choose grouping column (fall back to QC cluster if needed)
  clcol <- if ("cluster" %in% colnames(detected)) "cluster" else ".__qc_cluster__"
  if (!(clcol %in% colnames(detected))) {
    stop("rescue_by_coherence(): grouping cluster column not found.")
  }

  # ---- align rows across inputs
  rn  <- intersect(rownames(metrics), rownames(detected))
  det <- detected[rn, , drop = FALSE]
  met <- metrics[rn, , drop = FALSE]

  # ---- global medians for debug / gates (robust to missing columns)
  get_med <- function(v) if (is.null(v)) NA_real_ else stats::median(v, na.rm = TRUE)
  mt   <- if ("pctMT"        %in% colnames(met)) met$pctMT        else NULL
  nf   <- if ("nFeature"     %in% colnames(met)) met$nFeature     else NULL
  strs <- if ("stress_score" %in% colnames(met)) met$stress_score else NULL

  g_med <- list(
    mito   = get_med(mt),
    feat   = get_med(nf),
    stress = get_med(strs)
  )

  # ---- label coherence per cluster (dominant label fraction; Unknown=NA)
  labels <- obj@meta.data[rn, celltype_col, drop = TRUE]
  labels[is.na(labels) | labels == "Unknown"] <- NA
  cl <- det[[clcol]]

  tab      <- table(cl, labels, useNA = "no")
  tot      <- rowSums(tab)
  dom_frac <- ifelse(tot > 0, apply(tab, 1L, function(x) max(x) / sum(x)), 0)

  # ---- per-cluster medians and flags (robust to missing metrics)
  med_by <- function(v) {
    if (is.null(v)) return(rep(NA_real_, length(unique(cl))))
    tapply(v, cl, function(x) stats::median(x, na.rm = TRUE))
  }

  cluster_ids <- as.integer(names(tapply(cl, cl, length)))
  cm <- data.frame(
    cluster       = cluster_ids,
    size          = as.integer(tapply(cl, cl, length)),
    removed       = as.integer(tapply(det$qc_status == "remove", cl, sum)),
    med_mito      = as.numeric(med_by(mt)),
    med_feat      = as.numeric(med_by(nf)),
    med_stress    = as.numeric(med_by(strs)),
    stringsAsFactors = FALSE
  )
  rownames(cm) <- cm$cluster
  # handle coherence vector align
  cm$coherence <- as.numeric(dom_frac[as.character(cm$cluster)])
  cm$coherence[is.na(cm$coherence)] <- 0

  # gates (for debug table only; rescue is score-based below)
  cm$ok_size   <- cm$size >= as.integer(min_cluster_size)
  cm$pass_mito <- if (is.finite(g_med$mito))  cm$med_mito <= g_med$mito else TRUE
  cm$pass_feat <- if (is.finite(g_med$feat))  cm$med_feat >= g_med$feat else TRUE
  cm$pass_strs <- if (is.finite(g_med$stress)) cm$med_stress <= g_med$stress else TRUE

  # number of available metric gates and pass fraction
  cm$n_gates      <- as.integer((!is.na(cm$med_mito)) + (!is.na(cm$med_feat)) + (!is.na(cm$med_stress)))
  cm$n_pass       <- as.integer(cm$pass_mito) + as.integer(cm$pass_feat) + as.integer(cm$pass_strs)
  cm$pass_frac    <- ifelse(cm$n_gates > 0, cm$n_pass / cm$n_gates, 1.0)

  # ---------------------- continuous rescue score ---------------------------
  cm$removed_frac <- ifelse(cm$size > 0, cm$removed / cm$size, 1.0)
  cm$score        <- rowMeans(cbind(
    cm$coherence,         # 0..1, higher = better
    cm$pass_frac,         # 0..1, higher = better
    1 - cm$removed_frac   # 0..1, higher = better
  ), na.rm = TRUE)

  cm$pass_coh  <- cm$coherence >= coherence_min     # debug column
  cm$eligible  <- cm$ok_size & (cm$score >= resc_threshold)

  # rescue flags: only consider cells currently "remove" in eligible clusters
  rescue_flag <- rep(FALSE, nrow(det))
  eligible_clusters <- cm$cluster[cm$eligible %in% TRUE]
  if (length(eligible_clusters)) {
    ix <- which(det$qc_status == "remove" & (det[[clcol]] %in% eligible_clusters))
    if (length(ix)) rescue_flag[ix] <- TRUE
  }

  # reason strings (include score for transparency)
  reason <- character(nrow(det))
  reason[det$qc_status == "remove" & !rescue_flag] <- "not_eligible"
  if (length(eligible_clusters)) {
    scmap <- cm$score; names(scmap) <- cm$cluster
    rcl   <- det[[clcol]]
    reason[rescue_flag] <- sprintf("eligible_cluster(score=%.2f)", scmap[as.character(rcl[rescue_flag])])
  }

  # return frame aligned to 'det'
  out <- data.frame(
    rescue_flag   = rescue_flag[rownames(det)],
    rescue_reason = reason[rownames(det)],
    row.names     = rownames(det),
    check.names   = FALSE
  )

  # attach a rich gate table for external logging
  gate_tbl <- cm[order(-cm$removed), c(
    "cluster","size","removed","removed_frac","med_mito","med_feat","med_stress",
    "coherence","ok_size","pass_mito","pass_feat","pass_strs","n_pass","n_gates",
    "pass_frac","score","pass_coh","eligible"
  )]
  attr(out, "rescue_gate_table") <- gate_tbl

  # console trace
  if (getOption("scQCenrich.debug", FALSE)) {
    message(sprintf(
      "[rescue] mode=%s | resc_threshold=%.2f | min_cluster_size=%d | legacy(coh>=%.2f, gates_min_pass=%d)",
      rescue_mode, resc_threshold, as.integer(min_cluster_size), coherence_min, as.integer(gates_min_pass)
    ))
    message(sprintf(
      "[rescue] summary by cluster (top 5 by removed): %s",
      paste(utils::capture.output(print(utils::head(gate_tbl, 5))), collapse=" | ")
    ))
  }

  out
}
