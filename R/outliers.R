# ---- R/outliers.R ----
# Safe PCA embed + outlier detectors (GMM/HDBSCAN). No dense coercions, no zero-variance crashes.

# Build a robust QC feature matrix (impute, drop constants, robust-scale)
.qc_feature_matrix <- function(metrics) {
  stopifnot(is.data.frame(metrics))
  feats <- intersect(c("pctMT","nFeature","stress_score","MALAT1_frac",
                       "intronic_frac","u2s_ratio","spliced_frac","unspliced_frac"),
                     colnames(metrics))
  if (!length(feats)) stop("qc_embed(): no QC features found in metrics.")
  X <- as.data.frame(metrics[, feats, drop = FALSE])

  # Impute NAs with median (or 0 if entire column is NA)
  for (j in seq_along(feats)) {
    v <- X[[j]]; v[!is.finite(v)] <- NA
    if (all(is.na(v))) v[] <- 0 else v[is.na(v)] <- stats::median(v, na.rm = TRUE)
    X[[j]] <- v
  }

  # Drop constant / zero-variance columns (the root cause of prcomp error)
  uniq <- vapply(X, function(v) length(unique(v)), integer(1))
  sds  <- vapply(X, function(v) stats::sd(v), numeric(1))
  keep <- names(uniq)[uniq > 1 & is.finite(sds) & sds > 0]
  if (length(keep) < 2L) {
    # keep at least 1 column so the pipeline can proceed
    if (!length(keep)) keep <- feats[1]
  }
  X <- as.matrix(X[, keep, drop = FALSE])

  # Robust scale by median/MAD (avoid prcomp's internal rescale)
  med <- apply(X, 2, stats::median)
  mad <- apply(X, 2, stats::mad, constant = 1)
  mad[!is.finite(mad) | mad == 0] <- 1
  Z <- sweep(sweep(X, 2, med, "-"), 2, mad, "/")
  Z
}

# Low-d embedding of QC features (safe PCA)
qc_embed <- function(metrics, n_pcs = 10) {
  Z <- .qc_feature_matrix(metrics)
  if (ncol(Z) < 2L) {
    return(cbind(PC1 = as.numeric(scale(Z[,1]))))
  }
  pc <- stats::prcomp(Z, center = FALSE, scale. = FALSE, retx = TRUE)
  k <- min(ncol(pc$x), n_pcs)
  pc$x[, 1:k, drop = FALSE]
}

detect_outliers <- function(metrics,
                            method = c("gmm","hdbscan","threshold"),
                            sample_col = NULL,
                            remove_quantile = 0.90,
                            border_z = 1.5) {
  method <- match.arg(method)
  rn <- rownames(metrics); stopifnot(!is.null(rn))

  # ===== helpers ============================================================
  zfix <- function(v) {
    v <- as.numeric(v)
    s <- stats::sd(v, na.rm = TRUE)
    if (!is.finite(s) || s == 0) return(rep(0, length(v)))
    (v - mean(v, na.rm = TRUE)) / s
  }
  sdz <- function(v) {
    v <- as.numeric(v)
    m <- stats::median(v, na.rm=TRUE)
    s <- 1.4826 * stats::mad(v, constant=1, na.rm=TRUE)
    if (!is.finite(s) || s == 0) return(rep(0, length(v)))
    (v - m) / s
  }
  rZ   <- function(x) {  # robust z
    m <- stats::median(x, na.rm=TRUE)
    s <- 1.4826 * stats::mad(x, constant=1, na.rm=TRUE)
    if (!is.finite(s) || s == 0) return(rep(0, length(x)))
    (x - m) / s
  }
  rZabs <- function(x) abs(rZ(x))  # two-tailed deviation (penalize both high & low)

  get_pct <- function(x) {
    if (is.null(x)) return(NULL)
    if (all(is.na(x))) return(x)
    if (max(x, na.rm = TRUE) > 1) x <- x / 100
    x
  }

  # ---- extreme pctMT detector (unchanged logic) ---------------------------
  flag_extreme_pctmt <- function(df, pct, nfeat, sample) {
    message("detecting extreme pctmt")
    if (is.null(pct)) return(list(flag = rep(FALSE, nrow(df)), reason = rep("", nrow(df))))
    hard_cap <- getOption("scQCenrich.extreme_mito_cutoff",    0.50)
    pctl     <- getOption("scQCenrich.extreme_mito_percentile",0.995)
    nmads    <- getOption("scQCenrich.extreme_mito_nmads",     5)
    floor_p  <- getOption("scQCenrich.extreme_mito_floor",     0.50)

    lowN_feat <- rep(FALSE, nrow(df))
    if (!is.null(nfeat)) {
      q10 <- stats::quantile(nfeat, 0.10, na.rm = TRUE)
      lowN_feat <- nfeat <= q10
    }

    groups <- if (!is.null(sample) && sample %in% colnames(df)) df[[sample]] else factor("all")

    flag <- rep(FALSE, nrow(df)); reason <- rep("", nrow(df))
    for (g in levels(groups)) {
      idx <- which(groups == g); if (!length(idx)) next
      pg <- pct[idx]

      hard_hit <- pg >= hard_cap
      soft_thr <- stats::quantile(pg, probs = pctl, na.rm = TRUE)
      soft_hit <- pg >= soft_thr

      m  <- stats::median(pg, na.rm = TRUE)
      sd <- 1.4826 * stats::mad(pg, constant = 1, na.rm = TRUE)
      mad_thr <- max(floor_p, m + nmads * sd)
      mad_hit <- pg >= mad_thr

      flag[idx] <- hard_hit | ((soft_hit | mad_hit) & lowN_feat[idx])
      reason[idx] <- ifelse(hard_hit, "extreme_pctMT_hardcap",
                            ifelse((((soft_hit %in% TRUE) | (mad_hit %in% TRUE)) & (lowN_feat[idx] %in% TRUE)),
                                   "extreme_pctMT_supported_by_low_complexity",""))
    }
    list(flag = flag, reason = reason)
  }

  # ---- compact QC embedding (now includes intronic two-tailed dev) --------
  qc_embed <- function(df, n_pcs = 5) {
    cols <- c("pctMT","nFeature","stress_score","u2s_ratio")
    cols <- cols[cols %in% colnames(df)]
    if (length(cols) == 0L) stop("metrics lacks required QC columns.")

    X <- do.call(cbind, lapply(cols, function(nm) {
      x <- df[[nm]]
      if (nm == "nFeature") x <- -x
      zfix(x)
    }))

    # Add two-tailed intronic deviation as a dedicated axis (strong signal)
    if ("intronic_frac" %in% colnames(df)) {
      intr_dev <- rZabs(df$intronic_frac)
      X <- cbind(X, intr_dev)
      colnames(X)[ncol(X)] <- "intr_dev"
    }
    if ("MALAT1_frac" %in% colnames(df)) X <- cbind(X, zfix(df$MALAT1_frac))

    kp <- min(n_pcs, ncol(X))
    if (kp <= 1L) return(matrix(X[, 1L], ncol = 1L, dimnames = list(rownames(df), "PC1")))
    pcs <- stats::prcomp(X, center = FALSE, scale. = FALSE)$x
    pcs[, seq_len(min(kp, ncol(pcs))), drop = FALSE]
  }

  # ---- threshold mode (unchanged) -----------------------------------------
  if (method == "threshold") {
    mad_thr <- function(x) { m <- stats::median(x, na.rm = TRUE); s <- 1.4826 * stats::mad(x, constant = 1, na.rm = TRUE); list(low = m - 3*s, high = m + 3*s) }
    thr_mito   <- mad_thr(metrics$pctMT)
    thr_nfeat  <- mad_thr(metrics$nFeature)
    thr_stress <- mad_thr(metrics$stress_score)

    f_mito   <- metrics$pctMT        > thr_mito$high
    f_lowg   <- metrics$nFeature     < thr_nfeat$low
    f_stress <- metrics$stress_score > thr_stress$high

    status <- ifelse(f_mito | f_lowg | f_stress, "borderline", "keep")

    pfrac <- get_pct(metrics$pctMT)
    ext <- flag_extreme_pctmt(metrics, pfrac, metrics$nFeature, sample_col)
    ext_idx <- which(!is.na(ext$flag) & ext$flag)
    if (length(ext_idx)) status[ext_idx] <- "remove"

    reasons <- rep("", nrow(metrics))
    if (length(ext_idx)) {
      reasons[ext_idx] <- ext$reason[ext_idx]; reasons[is.na(reasons)] <- ""
    }
    return(data.frame(
      qc_status = factor(status, levels = c("keep","borderline","remove")),
      qc_score  = as.numeric(zfix(metrics$pctMT) + zfix(-metrics$nFeature)),
      cluster   = NA_integer_,
      reason    = reasons,
      row.names = rn, check.names = FALSE
    ))
  }

  # ---- model-based path (GMM/HDBSCAN) with intronic heavy two-tail --------
  emb <- qc_embed(metrics)

  used_fallback <- FALSE
  cl <- switch(method,
               gmm = {
                 if (!requireNamespace("mclust", quietly = TRUE)) stop("Package 'mclust' is required for method='gmm'.")
                 mns <- asNamespace("mclust")
                 out <- tryCatch({ mns$Mclust(emb, G = 1:6) }, error = function(e) {
                   if (grepl("mclustBIC", conditionMessage(e), ignore.case = TRUE)) { used_fallback <<- TRUE; NULL } else stop(e)
                 })
                 if (is.null(out)) {
                   if (!requireNamespace("dbscan", quietly = TRUE)) stop("mclust failed and 'dbscan' not installed for fallback.")
                   k <- max(20L, min(50L, round(0.01 * nrow(emb))))
                   as.integer(dbscan::hdbscan(emb, minPts = k)$cluster)
                 } else as.integer(out$classification)
               },
               hdbscan = {
                 if (!requireNamespace("dbscan", quietly = TRUE)) stop("Package 'dbscan' is required for method='hdbscan'.")
                 n <- nrow(emb); k <- max(20L, min(50L, round(0.01 * n)))
                 as.integer(dbscan::hdbscan(emb, minPts = k)$cluster)
               }
  )

  # cluster-level scoring (two-tailed intronic deviation, heavier weight)
  df <- metrics; df$.__cl__ <- cl

  # cell-level intronic two-tailed dev for later
  intr_dev_cell <- if ("intronic_frac" %in% colnames(metrics)) rZabs(metrics$intronic_frac) else rep(0, nrow(metrics))

  agg <- stats::aggregate(list(
    pctMT   = df$pctMT,
    negGene = -df$nFeature,
    stress  = df$stress_score,
    intrDev = intr_dev_cell,
    u2s     = df$u2s_ratio,
    malat   = if ("MALAT1_frac" %in% colnames(df)) df$MALAT1_frac else NA_real_
  ), by = list(cluster = df$.__cl__), FUN = function(v) mean(v, na.rm = TRUE))

  w_intr   <- getOption("scQCenrich.intr_two_tail_weight", 1.5)   # <- heavier default
  w_u2s    <- getOption("scQCenrich.u2s_weight", 1.0)

  zdf <- data.frame(
    mito = sdz(agg$pctMT),
    gene = sdz(agg$negGene),
    strs = sdz(agg$stress),
    intr = w_intr * sdz(agg$intrDev),   # two-tailed intronic deviation, up-weighted
    u2s  = w_u2s  * sdz(agg$u2s),
    mal  = if (all(is.na(agg$malat))) 0 else 0.5 * sdz(agg$malat)
  )
  cl_score <- rowSums(zdf, na.rm = TRUE)

  thr   <- stats::quantile(cl_score, probs = remove_quantile, na.rm = TRUE)
  worst <- unique(c(agg$cluster[cl_score >= thr], agg$cluster[which.max(cl_score)]))

  status <- ifelse(cl %in% worst, "remove", "keep")

  # borderline rules (include intronic two-tailed dev)
  mito_z <- zfix(df$pctMT)
  gene_z <- zfix(-df$nFeature)
  intr_border <- getOption("scQCenrich.intr_border_z", 2)
  intr_escal  <- getOption("scQCenrich.intr_escalate_z", 3)

  borderline <- ((mito_z > border_z) | (gene_z > border_z)) & !(cl %in% worst)
  # add intronic two-tailed deviation condition
  borderline <- borderline | (intr_dev_cell > intr_border & !(cl %in% worst))

  # extreme pctMT override & escalation
  pfrac <- get_pct(metrics$pctMT)
  ext <- flag_extreme_pctmt(metrics, pfrac, metrics$nFeature, sample_col)
  ext_idx <- which(!is.na(ext$flag) & ext$flag)
  if (length(ext_idx)) status[ext_idx] <- "remove"

  esc_cap <- getOption("scQCenrich.escalate_mito_cutoff", 0.40)
  esc_z   <- getOption("scQCenrich.escalate_mito_z", 3.0)
  esc_idx <- which(status == "borderline" &
                     !is.na(pfrac) & (pfrac >= esc_cap | mito_z >= esc_z))
  if (length(esc_idx)) status[esc_idx] <- "remove"

  # intronic escalation: very extreme two-tailed intronic deviation -> remove
  intr_esc_idx <- which(status == "borderline" & intr_dev_cell >= intr_escal)
  if (length(intr_esc_idx)) status[intr_esc_idx] <- "remove"

  # ---- SAFETY VALVE: cap overall removal fraction (unchanged) -------------
  cap <- getOption("scQCenrich.max_remove_fraction", 0.1)
  changed_idx <- integer(0)

  if (mean(status == "remove", na.rm = TRUE) > cap) {
    used_safety <- FALSE
    qs <- seq(min(remove_quantile + 0.02, 0.99), 0.99, by = 0.02)
    for (q in qs) {
      thr2   <- stats::quantile(cl_score, probs = q, na.rm = TRUE)
      worst2 <- unique(c(agg$cluster[cl_score >= thr2], agg$cluster[which.max(cl_score)]))
      status2 <- ifelse(cl %in% worst2, "remove", "keep")
      borderline2 <- ((mito_z > border_z) | (gene_z > border_z)) & !(cl %in% worst2)
      borderline2 <- borderline2 | (intr_dev_cell > intr_border & !(cl %in% worst2))

      if (length(ext_idx)) status2[ext_idx] <- "remove"
      esc_idx2 <- which(status2 == "borderline" &
                          !is.na(pfrac) & (pfrac >= esc_cap | mito_z >= esc_z))
      if (length(esc_idx2)) status2[esc_idx2] <- "remove"
      intr_esc_idx2 <- which(status2 == "borderline" & intr_dev_cell >= intr_escal)
      if (length(intr_esc_idx2)) status2[intr_esc_idx2] <- "remove"

      if (mean(status2 == "remove", na.rm = TRUE) <= cap) {
        status <- status2
        used_safety <- TRUE
        break
      }
    }

    if (!used_safety && mean(status == "remove", na.rm = TRUE) > cap) {
      rem_idx <- which(status == "remove")
      if (length(rem_idx)) {
        sev <- pmax(mito_z[rem_idx], 0) + pmax(gene_z[rem_idx], 0) +  # tie-breaker includes intronic dev
          pmax(intr_dev_cell[rem_idx], 0)
        target_n <- ceiling(length(status) * cap)
        n_demote <- max(0L, length(rem_idx) - target_n)
        if (n_demote > 0L) {
          demote <- rem_idx[order(sev, decreasing = FALSE)[seq_len(n_demote)]]
          status[demote] <- "borderline"
          changed_idx <- demote
        }
      }
    }
  }

  # ---- reasons and return --------------------------------------------------
  reasons <- rep("", nrow(df))
  if (length(ext_idx)) {
    reasons[ext_idx] <- ext$reason[ext_idx]
    reasons[is.na(reasons)] <- ""
  }
  if (length(esc_idx)) {
    reasons[esc_idx] <- ifelse(nzchar(reasons[esc_idx]),
                               paste0(reasons[esc_idx], "; escalated_pctMT"),
                               "escalated_pctMT")
  }
  if (length(intr_esc_idx)) {
    add <- "escalated_intronic_two_tail"
    reasons[intr_esc_idx] <- ifelse(nzchar(reasons[intr_esc_idx]),
                                    paste0(reasons[intr_esc_idx], "; ", add), add)
  }
  if (length(changed_idx)) {
    add <- "safety_cap_lenient_throttle"
    reasons[changed_idx] <- ifelse(nzchar(reasons[changed_idx]),
                                   paste0(reasons[changed_idx], "; ", add), add)
  }

  if (getOption("scQCenrich.debug", FALSE)) {
    msg <- sprintf("[detect_outliers] worst=%s | remove=%d borderline=%d keep=%d | thr=%.3f border_z=%.2f | w_intr=%.2f intr_border=%.2f intr_escal=%.2f",
                   paste(sort(unique(worst)), collapse=","), sum(status=="remove"),
                   sum(status=="borderline"), sum(status=="keep"), as.numeric(thr),
                   border_z, w_intr, intr_border, intr_escal)
    if (used_fallback) msg <- paste0(msg, " | NOTE: GMM -> HDBSCAN fallback (mclustBIC)")
    hard_n <- sum(get_pct(metrics$pctMT) >= getOption("scQCenrich.extreme_mito_cutoff", 0.50), na.rm = TRUE)
    if (hard_n > 0) msg <- paste0(msg, sprintf(" | extreme_pctMT_hardcap=%d", hard_n))
    message(msg)
  }

  data.frame(
    qc_status = factor(ifelse(borderline & status=="keep","borderline", status),
                       levels = c("keep","borderline","remove")),
    qc_score  = as.numeric(zfix(metrics$pctMT) + zfix(-metrics$nFeature) + rZabs(metrics$intronic_frac)),
    cluster   = cl,
    reason    = ifelse(borderline & status=="keep", "high_mito_or_low_genes_or_intronic_dev",
                       ifelse(status=="remove",
                              ifelse(nzchar(reasons), reasons, "worst_cluster_by_QC"),
                              "")),
    row.names = rn, check.names = FALSE
  )
}
