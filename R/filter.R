#' Flag low-quality cells
#' @param metrics data.frame of QC metrics (rows=cells). If it contains a logical
#'   column 'is_doublet' (as added by calcQCmetrics), it will be honored.
#' @param method "threshold", "gmm", or "hdbscan"
#' @param score_cutoff integer for threshold mode (default 2)
#' @param auto_k kept for API compat (unused here)
#' @param obj optional Seurat object for rescue step
#' @param celltype_col optional meta column for rescue label coherence
#' @param sample_col optional sample column (passed to detector; may be NULL)
#' @param rescue_mode "moderate","conservative","aggressive"
#' @param min_cluster_size minimum cluster size to consider for rescue; NULL = adaptive (~0.5% or >=10)
#' @param doublet_action what to do with predicted doublets in `metrics$is_doublet`:
#'   "remove" (default), "borderline", or "none"
#' @param doublet_col name of the doublet flag column in `metrics` (default "is_doublet")
#' @param qc_strength One of c("auto","default","lenient","strict")
#' @param save_dir Character path or NULL. If non-NULL, write intermediate CSVs here.
#' @return data.frame with rownames=cells and columns: qc_status, qc_score, cluster, reason,
#' @export
flagLowQuality <- function(
    metrics,
    method = c("threshold","gmm","hdbscan"),
    score_cutoff = 2,
    auto_k = 3,
    obj = NULL,
    celltype_col = NULL,
    sample_col = NULL,
    rescue_mode = c("moderate","conservative","aggressive"),
    min_cluster_size = NULL,   # NULL -> adaptive
    doublet_action = c("remove","borderline","none"),
    doublet_col = "is_doublet",
    qc_strength = c("auto","default","lenient","strict"),
    save_dir = NULL
){
  qc_strength    <- match.arg(qc_strength)
  method         <- match.arg(method)
  rescue_mode    <- match.arg(rescue_mode)
  doublet_action <- match.arg(doublet_action)

  s <- switch(qc_strength,
              lenient = list(remove_quantile = 0.98, border_z = 2.0, threshold_score = 3L, rescue = "aggressive"),
              strict  = list(remove_quantile = 0.85, border_z = 0.8, threshold_score = 1L, rescue = "conservative"),
              default = list(remove_quantile = 0.95, border_z = 1.2, threshold_score = 2L, rescue = "moderate"),
              auto    = list(remove_quantile = NA_real_, border_z = 1.2, threshold_score = 2L, rescue = "moderate")
  )

  ## --- Learn remove_quantile in 'auto' via a quick pre-score
  if (identical(qc_strength, "auto")) {
    z <- function(x, sign = +1L) {
      if (is.null(x) || all(!is.finite(x))) return(rep(0, nrow(metrics)))
      zx <- as.numeric(scale(x)); zx[!is.finite(zx)] <- 0
      sign * zx
    }
    pre <- z(metrics$pctMT, +1L) + z(metrics$nFeature, -1L) + z(metrics$stress_score, +1L)
    if ("intronic_frac" %in% colnames(metrics)) pre <- pre + z(metrics$intronic_frac, +1L)

    rq <- tryCatch({
      qs <- stats::quantile(pre, probs = c(0.65, 0.75, 0.85), na.rm = TRUE, type = 7)
      thr <- qs[[2L]]
      if (is.finite(thr)) stats::ecdf(pre)(thr) else NA_real_
    }, error = function(e) NA_real_)
    if (!is.finite(rq)) rq <- 0.95
    rq <- max(0.85, min(0.97, rq))
    s$remove_quantile <- rq
    if (isTRUE(getOption("scQCenrich.debug", FALSE))) {
      message(sprintf("[qc_strength=auto] learned remove_quantile=%.3f (target removal ~%.0f%%)",
                      s$remove_quantile, (1 - s$remove_quantile)*100))
    }
  }

  .robust_limits <- function(x, nmads = 3, invert = FALSE) {
    if (is.null(x) || all(!is.finite(x))) return(list(low=NA_real_, high=NA_real_))
    xx <- as.numeric(x)
    m  <- stats::median(xx, na.rm = TRUE)
    s  <- stats::mad(xx, constant = 1.4826, na.rm = TRUE)
    if (!invert) list(low = m - nmads * s, high = m + nmads * s)
    else         list(low = m + nmads * s, high = m - nmads * s)
  }

  .apply_doublet <- function(df) {
    if (!doublet_col %in% colnames(df)) return(df)
    isdb <- df[[doublet_col]]
    if (!is.logical(isdb)) isdb <- as.logical(isdb)
    if (!any(isTRUE(isdb))) return(df)
    idx <- which(isTRUE(isdb))
    if (doublet_action == "remove") {
      df$qc_status[idx] <- factor("remove", levels = c("keep","borderline","remove"))
      add <- "doublet_detected"
      df$reason[idx] <- ifelse(nzchar(df$reason[idx]),
                               paste0(df$reason[idx], "; ", add), add)
    } else if (doublet_action == "borderline") {
      df$qc_status[idx] <- factor("borderline", levels = c("keep","borderline","remove"))
      add <- "doublet_detected"
      df$reason[idx] <- ifelse(nzchar(df$reason[idx]),
                               paste0(df$reason[idx], "; ", add), add)
    }
    df
  }

  det <- detect_outliers(
    metrics         = metrics,
    method          = if (method == "threshold") "threshold" else if (method=="hdbscan") "hdbscan" else "gmm",
    sample_col      = sample_col,
    remove_quantile = s$remove_quantile,
    border_z        = s$border_z
  )
  # ensure rownames present and aligned
  if (is.null(rownames(det))) rownames(det) <- rownames(metrics)

  .auto_min_scater_frac <- 0.80
  .auto_scater_enforce  <- TRUE
  if (identical(qc_strength, "auto") && isTRUE(.auto_scater_enforce)) {

    .scater_like_drop <- function(df, nmads = 3L) {
      pick1 <- function(cands) { cands[cands %in% names(df)][1] }
      lib_col  <- pick1(c("nCount","nCount_RNA","sum","libsize","total_counts","UMI"))
      feat_col <- pick1(c("nFeature","nFeature_RNA","detected","n_genes","nGenes"))
      mito_col <- pick1(c("pctMT","percent.mt","mito_frac","mt_frac","subsets_Mt_percent"))
      n <- nrow(df)
      if (is.na(lib_col) || is.na(feat_col) || is.na(mito_col) || n == 0L) return(rep(FALSE, n))
      lib  <- as.numeric(df[[lib_col]])
      feat <- as.numeric(df[[feat_col]])
      mito <- as.numeric(df[[mito_col]])
      if (all(is.finite(mito)) && any(mito > 1, na.rm = TRUE)) mito <- mito / 100
      thr_lib   <- .robust_limits(lib,  nmads = nmads)
      thr_feat  <- .robust_limits(feat, nmads = nmads)
      thr_mito  <- .robust_limits(mito, nmads = nmads)
      low_lib   <- lib  < thr_lib$low
      low_feat  <- feat < thr_feat$low
      high_mito <- mito > thr_mito$high
      (low_lib | low_feat | high_mito)    # vector; no isTRUE() on a vector
    }

    scater_drop <- .scater_like_drop(metrics, nmads = 3L)
    n_scater    <- sum(scater_drop, na.rm = TRUE)
    target_n    <- ceiling(.auto_min_scater_frac * n_scater)
    cur_remove_n<- sum(det$qc_status == "remove", na.rm = TRUE)

    if (is.finite(target_n) && target_n > 0L && cur_remove_n < target_n) {
      z <- function(x, sign = +1L) {
        if (is.null(x) || all(!is.finite(x))) return(rep(0, nrow(metrics)))
        zx <- as.numeric(scale(x)); zx[!is.finite(zx)] <- 0
        sign * zx
      }
      pre <- z(metrics$pctMT, +1L) + z(metrics$nFeature, -1L) + z(metrics$stress_score, +1L)
      if ("intronic_frac" %in% colnames(metrics)) pre <- pre + z(metrics$intronic_frac, +1L)

      cand <- which(det$qc_status != "remove")
      need <- target_n - cur_remove_n
      if (length(cand) && need > 0L) {
        cand <- cand[order(pre[cand], decreasing = TRUE)]
        take <- head(cand, need)
        det$qc_status[take] <- factor("remove", levels = c("keep","borderline","remove"))
        add <- sprintf("auto_min_scater_enforce(>=%.0f%% of scater)", .auto_min_scater_frac*100)
        det$reason[take] <- ifelse(nzchar(det$reason[take]),
                                   paste0(det$reason[take], "; ", add), add)
      }
    }
  }

  if (identical(method, "threshold")) {
    thr_mito   <- .robust_limits(metrics$pctMT,        nmads = 3L)
    thr_nfeat  <- .robust_limits(metrics$nFeature,     nmads = 3L)
    thr_stress <- .robust_limits(metrics$stress_score, nmads = 3L)
    pmt <- metrics$pctMT
    if (all(is.finite(pmt)) && any(pmt > 1, na.rm = TRUE)) pmt <- pmt / 100
    f1 <- ifelse(is.na(pmt),                  FALSE, pmt                 > thr_mito$high)
    f2 <- ifelse(is.na(metrics$nFeature),     FALSE, metrics$nFeature     < thr_nfeat$low)
    f5 <- ifelse(is.na(metrics$stress_score), FALSE, metrics$stress_score > thr_stress$high)

    score  <- as.integer(f1) + as.integer(f2) + as.integer(f5)
    status <- ifelse(score >= as.integer(s$threshold_score), "remove",
                     ifelse(score > 0L, "borderline", "keep"))

    det$qc_status <- factor(status, levels = c("keep","borderline","remove"))
    det$qc_score  <- as.numeric(scale(metrics$pctMT, center = TRUE, scale = TRUE)) +
      as.numeric(scale(-metrics$nFeature, center = TRUE, scale = TRUE)) +
      if ("intronic_frac" %in% names(metrics)) abs(as.numeric(scale(metrics$intronic_frac)))
    det$reason    <- ifelse(det$qc_status=="remove","threshold_rule",
                            ifelse(det$qc_status=="borderline","one_metric_flag",""))
  }

  if ("pctMT" %in% colnames(metrics) && "intronic_frac" %in% colnames(metrics)) {
    rn   <- rownames(metrics)
    pmt  <- as.numeric(metrics[rn, "pctMT", drop = TRUE])
    intr <- as.numeric(metrics[rn, "intronic_frac", drop = TRUE])
    if (any(is.finite(pmt)) && any(is.finite(intr))) {
      if (max(pmt, na.rm = TRUE) > 1) pmt <- pmt / 100
      nuclei_cut <- getOption("scQCenrich.nuclei_median_intronic", 0.40)  # nuclei if median intronic ≥ 0.40
      med_intr   <- suppressWarnings(stats::median(intr, na.rm = TRUE))
      is_nuclei  <- is.finite(med_intr) && med_intr >= nuclei_cut

      mito_susp_cut   <- getOption("scQCenrich.mito_suspicious_cut", 0.25)  # ≥25% suspicious
      mito_high_cut   <- getOption("scQCenrich.mito_high_cut",       0.35)  # ≥35% remove
      intr_normal_max <- getOption("scQCenrich.normal_intronic_max", 0.35)  # ≤35% = normal/mid (scRNA)

      if (!is_nuclei) {
        idx_guard  <- which(is.finite(pmt) & is.finite(intr) & pmt >= mito_susp_cut & intr <= intr_normal_max)
        if (length(idx_guard)) {
          idx_remove <- idx_guard[pmt[idx_guard] >= mito_high_cut]
          idx_border <- setdiff(idx_guard, idx_remove)

          if (length(idx_border)) {
            set_to_border <- idx_border[det$qc_status[idx_border] == "keep"]
            if (length(set_to_border)) {
              det$qc_status[set_to_border] <- factor("borderline", levels = c("keep","borderline","remove"))
              add <- "mito_norm_intronic_guard"
              det$reason[set_to_border] <- ifelse(nzchar(det$reason[set_to_border]),
                                                  paste0(det$reason[set_to_border], "; ", add), add)
            }
          }
          if (length(idx_remove)) {
            to_rm <- idx_remove[det$qc_status[idx_remove] != "remove"]
            if (length(to_rm)) {
              det$qc_status[to_rm] <- factor("remove", levels = c("keep","borderline","remove"))
              add <- "mito_norm_intronic_guard_strict"
              det$reason[to_rm] <- ifelse(nzchar(det$reason[to_rm]),
                                          paste0(det$reason[to_rm], "; ", add), add)
            }
          }
        }
      }
    }
  }

  if (!"rescue_flag" %in% names(det))     det$rescue_flag    <- FALSE
  if (!"border_suggest" %in% names(det))  det$border_suggest <- FALSE

  if (!is.null(obj) && inherits(obj, "Seurat") && !is.null(celltype_col)) {
    if (!isTRUE(celltype_col %in% colnames(obj@meta.data))) {
      stop("flagLowQuality(): `celltype_col` not found in obj@meta.data.")
    }
    # ensure QC fields exist in meta for rescue gating (pull from metrics if missing)
    md <- obj@meta.data
    keymap <- list(
      pctMT_QC         = "pctMT",
      nFeature_QC      = "nFeature",
      stress_score_QC  = "stress_score",
      intronic_frac_QC = "intronic_frac",
      u2s_ratio_QC     = "u2s_ratio"
    )
    for (nm in names(keymap)) {
      src <- keymap[[nm]]
      if (!(nm %in% colnames(md)) && (src %in% colnames(metrics))) {
        mm <- metrics[rownames(md), src, drop = TRUE]
        md[, nm] <- as.numeric(mm)
      }
    }
    obj@meta.data <- md

    mcs <- if (!is.null(min_cluster_size) && is.finite(min_cluster_size)) as.integer(min_cluster_size) else 50L

    resc <- tryCatch(
      rescue_by_coherence(
        obj, metrics, det,
        celltype_col     = celltype_col,
        min_cluster_size = mcs,
        rescue_mode      = rescue_mode
      ),
      error = function(e) NULL
    )

    if (!is.null(resc) && NROW(resc) > 0L) {
      if (is.null(rownames(resc))) {
        id_col <- c("cell","barcode","Cell")[c("cell","barcode","Cell") %in% names(resc)][1]
        if (!is.na(id_col)) rownames(resc) <- as.character(resc[[id_col]])
      }
      if (is.null(rownames(resc))) {
        warning("[rescue] table has no rownames or cell column; skipping aggregation.")
      } else {
        g  <- rownames(resc)
        rf <- if ("rescue_flag"    %in% names(resc)) tapply(as.logical(resc[["rescue_flag"]]),    g, function(v) any(isTRUE(v))) else NULL
        bs <- if ("border_suggest" %in% names(resc)) tapply(as.logical(resc[["border_suggest"]]), g, function(v) any(isTRUE(v))) else NULL

        if (!is.null(rf)) {
          m <- match(names(rf), rownames(det)); ok <- which(!is.na(m))
          if (length(ok)) det$rescue_flag[m[ok]] <- as.logical(rf[ok])
        }
        if (!is.null(bs)) {
          m <- match(names(bs), rownames(det)); ok <- which(!is.na(m))
          if (length(ok)) det$border_suggest[m[ok]] <- as.logical(bs[ok])
        }

        if (isTRUE(getOption("scQCenrich.debug", FALSE))) {
          gt <- attr(resc, "rescue_gate_table")
          if (!is.null(gt)) {
            outdir <- if (!is.null(save_dir) && nzchar(save_dir)) save_dir else "qc_outputs"
            if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
            utils::write.csv(gt, file.path(outdir, "rescue_gate_table.csv"), row.names = FALSE)
            message(sprintf("[flagLowQuality] wrote %s (clusters=%d)",
                            file.path(outdir, "rescue_gate_table.csv"),
                            length(unique(gt$cluster_id))))
          }
        }
      }
    } else if (isTRUE(getOption("scQCenrich.debug", FALSE))) {
      message("[rescue] skipped (no eligible clusters or empty rescue table).")
    }


    ix <- which(det$rescue_flag %in% TRUE &
                  det$qc_status == "remove" &
                  !grepl("(?:^|;)\\s*(extreme_pctMT|mito_norm_intronic_guard_strict)\\b",
                         det$reason, perl = TRUE))
    if (length(ix)) {
      det$qc_status[ix] <- factor("borderline", levels = c("keep","borderline","remove"))
      det$reason[ix]    <- ifelse(nzchar(det$reason[ix]),
                                  paste0(det$reason[ix], "; rescued_by_coherence"),
                                  "rescued_by_coherence")
    }
  }

  det <- .apply_doublet(det)

  if ("pctMT" %in% colnames(metrics)) {
    pfrac <- metrics[rownames(det), "pctMT", drop = TRUE]
    if (max(pfrac, na.rm = TRUE) > 1) pfrac <- pfrac / 100
    cap <- getOption("scQCenrich.extreme_mito_cutoff", 0.50)
    viol <- which(!is.na(pfrac) & pfrac >= cap & det$qc_status != "remove")
    if (length(viol)) {
      det$qc_status[viol] <- factor("remove", levels = c("keep","borderline","remove"))
      add <- "enforced_extreme_pctMT_final"
      det$reason[viol] <- ifelse(nzchar(det$reason[viol]),
                                 paste0(det$reason[viol], "; ", add), add)
      if (isTRUE(getOption("scQCenrich.debug", FALSE))) {
        outdir <- if (!is.null(save_dir) && nzchar(save_dir)) save_dir else "qc_outputs"
        if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
        utils::write.csv(
          data.frame(barcode = rownames(det)[viol], pctMT = pfrac[viol], reason = det$reason[viol]),
          file.path(outdir, "extreme_mito_forced_removals.csv"), row.names = FALSE
        )
      }
    }
  }

  det
}
