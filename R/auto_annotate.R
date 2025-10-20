#' Auto-annotate cell types (multi-backend wrapper)
#'
#' @param obj Seurat object.
#' @param species Character, "mouse" or "human".
#' @param assay Assay to use (default "RNA").
#' @param prefer Optional character; preferred backend(s) to try first.
#' @param level "cell" or "cluster".
#' @param marker_score_level Level for marker_score ("cell"/"cluster").
#' @param custom_signatures Optional list/data.frame of custom markers.
#' @param ref_seurat Optional Seurat reference object (label transfer).
#' @param ref_label_col Column in `ref_seurat@meta.data` with labels.
#' @param unknown_min_genes Integer; minimum genes to avoid "unknown".
#' @param unknown_min_margin Numeric; min score margin to call a label.
#' @param unknown_min_top Integer; min top markers expressed to call a label.
#' @param organ Optional organ/tissue context string.
#' @param sctype_db_path Optional path to ScType flat DB (xlsx).
#' @param marker_method Marker backend when using marker-based annotation.
#' @param panglao_file Optional PanglaoDB CSV/XLSX file path.
#' @param panglao_tissue Optional tissue filter for PanglaoDB.
#' @param panglao_ui_max Integer; cap markers per type to avoid overplot.
#' @param panglao_canonical_only Logical; keep only canonical markers.
#' @param fm_top_n Integer; FindMarkers top N to keep per cluster.
#' @param fm_padj_max Numeric; padj cutoff for FindMarkers filtering.
#' @param overlap_mode How to combine/overlap multiple marker sources.
#' @return Seurat object with annotation columns added.
#' @export
auto_annotate <- function(obj,
                          species = c("mouse", "human"),
                          assay   = "RNA",
                          prefer  = c("sctype","marker_score","singleR","label_transfer"),
                          level   = c("cluster", "cell"),
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
                          panglao_file = system.file("extdata","PanglaoDB_markers_27_Mar_2020.tsv",
                                                     package = utils::packageName()),
                          panglao_tissue = NULL,
                          panglao_ui_max = 0.20,
                          panglao_canonical_only = TRUE,
                          fm_top_n = 50,
                          fm_padj_max = 0.05,
                          overlap_mode = c("count", "jaccard")) {
  stopifnot(inherits(obj, "Seurat"))
  species <- match.arg(species)
  prefer  <- match.arg(prefer, several.ok = TRUE)
  level   <- match.arg(level)
  marker_score_level <- match.arg(marker_score_level)
  marker_method <- match.arg(marker_method)
  overlap_mode  <- match.arg(overlap_mode)

  core_dir <- .aa_dbg_dir("core")
  .aa_log(core_dir,
          "[auto_annotate] start | species=%s assay=%s prefer=%s level=%s ms_level=%s method=%s",
          species, assay, paste(prefer, collapse="->"), level, marker_score_level, marker_method)
  .aa_log(core_dir, sprintf("cells=%d features=%d",
                            ncol(obj), nrow(Seurat::GetAssayData(obj, assay = assay, layer = "counts"))))

  backends <- unique(c(prefer, setdiff(c("singleR","label_transfer","marker_score","sctype"), prefer)))
  .dbg("auto_annotate(): backends (in order) = %s", paste(backends, collapse = " \u2192 "))

  # helper to write results -> meta.data + snippet + debug bundle
  .write_res <- function(obj, res, backend_name, method_tag = NULL) {
    stopifnot(is.list(res))
    dbgdir <- .aa_dbg_dir(backend_name)
    # cell labels & margins
    labs_vec <- res$labels
    conf_vec <- res$margin

    # ensure clusters if we need a cluster fallback
    if (anyNA(labs_vec)) {
      if (!"seurat_clusters" %in% colnames(obj@meta.data)) {
        obj <- .ensure_clusters(obj, assay)    # will create seurat_clusters if missing
      }
      if (!is.null(res$cluster_labels)) {
        cl <- as.character(obj$seurat_clusters); names(cl) <- colnames(obj)
        na_cells <- names(labs_vec)[is.na(labs_vec)]
        if (length(na_cells)) labs_vec[na_cells] <- res$cluster_labels[ cl[na_cells] ]
      }
      labs_vec[is.na(labs_vec)] <- "Unknown"
    }
    if (is.null(names(conf_vec)) && length(conf_vec) == length(labs_vec)) {
      names(conf_vec) <- names(labs_vec)
    }

    snip <- .annot_log_snippet(obj, list(labels = labs_vec,
                                         cluster_labels = res$cluster_labels,
                                         level = res$level),
                               backend = backend_name)
    cat(paste(snip, collapse = "\n"), "\n",
        file = file.path(dbgdir, "RESULT_SNIPPET.txt"), append = TRUE)
    .aa_log(dbgdir, paste(snip, collapse = "\n"))

    # write into meta.data
    src_val <- method_tag %||% backend_name
    cell_ids <- Seurat::Cells(obj)
    nm <- intersect(cell_ids, names(labs_vec))
    md <- obj@meta.data
    if (!("auto_celltype"      %in% colnames(md))) md$auto_celltype      <- NA_character_
    if (!("auto_celltype_src"  %in% colnames(md))) md$auto_celltype_src  <- NA_character_
    if (!("auto_celltype_conf" %in% colnames(md))) md$auto_celltype_conf <- NA_real_
    md[nm, "auto_celltype"]      <- as.character(labs_vec[nm])
    md[nm, "auto_celltype_src"]  <- src_val
    md[nm, "auto_celltype_conf"] <- suppressWarnings(as.numeric(conf_vec[nm]))
    obj@meta.data <- md
    .aa_log(dbgdir, "[auto_annotate] wrote auto_celltype for %d/%d cells",
            sum(!is.na(obj@meta.data[nm, "auto_celltype"])), length(nm))

    # keep debug bundle if present
    dbg <- attr(res, "debug")
    if (!is.null(dbg)) obj@misc$scQCenrich[[paste0(backend_name,"_debug")]] <- dbg
    obj
  }

  for (bk in backends) {
    if (bk == "sctype") {
      dir <- .aa_dbg_dir("sctype")
      .aa_log(dir, "[sctype] begin (tissues=%s, db=%s)",
              paste(organ %||% character(0), collapse="|"), sctype_db_path %||% "<none>")
      res <- try(.annot_sctype(obj, species, assay, level, tissues = organ, db_path = sctype_db_path), silent = FALSE)
      if (!inherits(res, "try-error") && !is.null(res)) {
        obj <- .write_res(obj, list(labels = res$cell_labels, margin = res$conf, cluster_labels = res$cluster_labels, level = level),
                          backend_name = "sctype")
        return(obj)
      }
      .aa_log(dir, "[sctype] no labels assigned; continuing")
    }

    if (bk == "singleR" && .has_singleR()) {
      dir <- .aa_dbg_dir("singler")
      .aa_log(dir, "[singleR] begin (level=%s)", level)
      res <- try(.annot_singleR(obj, species = species, assay = assay, level = level), silent = FALSE)
      if (!inherits(res, "try-error") && !is.null(res)) {
        obj <- .write_res(obj, list(labels = res$cell_labels, margin = res$conf, cluster_labels = NULL, level = level),
                          backend_name = "singleR")
        return(obj)
      }
      .aa_log(dir, "[singleR] failed or empty; continuing")
    }

    if (bk == "label_transfer" && !is.null(ref_seurat) && !is.null(ref_label_col)) {
      dir <- .aa_dbg_dir("label_transfer")
      .aa_log(dir, "[label_transfer] begin (ref_label_col=%s)", ref_label_col)
      res <- try(.annot_label_transfer(obj, ref_seurat, ref_label_col, assay = assay), silent = FALSE)
      if (!inherits(res, "try-error") && !is.null(res)) {
        labs <- res; names(labs) <- names(res)
        obj <- .write_res(obj, list(labels = labs, margin = setNames(rep(NA_real_, length(labs)), names(labs)),
                                    cluster_labels = NULL, level = level),
                          backend_name = "label_transfer")
        return(obj)
      }
      .aa_log(dir, "[label_transfer] failed or empty; continuing")
    }

    if (bk == "marker_score") {
      dir <- .aa_dbg_dir(if (identical(marker_method, "findmarkers")) "panglao_overlap" else "marker_score")
      .aa_log(dir, "[marker_score] begin (method=%s, level=%s)", marker_method, marker_score_level)
      res <- try(.annot_marker_score_generic(
        obj,
        species = species,
        assay = assay,
        level = marker_score_level,
        custom_signatures = custom_signatures,
        unknown_min_genes  = unknown_min_genes,
        unknown_min_margin = unknown_min_margin,
        unknown_min_top    = unknown_min_top,
        method = marker_method,
        panglao_file = panglao_file,
        panglao_tissue = panglao_tissue,
        panglao_ui_max = panglao_ui_max,
        panglao_canonical_only = panglao_canonical_only,
        fm_top_n = fm_top_n,
        fm_padj_max = fm_padj_max,
        overlap_mode = overlap_mode,
        debug = TRUE
      ), silent = FALSE)

      if (!inherits(res, "try-error") && !is.null(res)) {
        obj <- .write_res(
          obj,
          list(labels = res$labels, margin = res$margin, cluster_labels = res$cluster_labels, level = res$level),
          backend_name = if (identical(marker_method, "findmarkers")) "panglao_overlap" else "marker_score",
          method_tag   = if (identical(marker_method, "findmarkers")) "panglao_overlap" else "marker_score"
        )
        return(obj)
      }
      .aa_log(dir, "[marker_score] failed or empty; continuing")
    }
  }

  warning("auto_annotate(): no backend succeeded; no labels assigned.")
  obj
}



.has_singleR <- function() {
  requireNamespace("SingleR", quietly = TRUE) &&
    requireNamespace("celldex", quietly = TRUE) &&
    requireNamespace("SingleCellExperiment", quietly = TRUE)
}

.ensure_clusters <- function(obj, assay) {
  dir <- .aa_dbg_dir("core")
  Seurat::DefaultAssay(obj) <- assay

  # If clusters already there, reuse
  if ("seurat_clusters" %in% colnames(obj@meta.data)) {
    .aa_log(dir, "[prep] using existing seurat_clusters; K=%d",
            length(unique(obj$seurat_clusters)))
    return(obj)
  }

  .aa_log(dir, "[prep] computing clusters (assay=%s)", assay)

  # Ensure normalized 'data' + HVGs
  obj <- Seurat::NormalizeData(obj, verbose = FALSE)
  obj <- Seurat::FindVariableFeatures(obj, nfeatures = 2000, verbose = FALSE)
  hvgs <- Seurat::VariableFeatures(obj)

  # Scale **only HVGs** (explicitly)
  if (length(hvgs)) {
    obj <- Seurat::ScaleData(obj, features = hvgs, verbose = FALSE)
  } else {
    # degenerate case; fall back to all features to avoid empty PCA
    obj <- Seurat::ScaleData(obj, verbose = FALSE)
  }

  # PCA on HVGs; cap npcs by #features available
  npcs <- min(30L, max(2L, length(hvgs)))
  obj  <- Seurat::RunPCA(obj, features = hvgs, npcs = npcs, verbose = FALSE)

  # Dims for neighbors must not exceed available PCs
  nd <- min(20L, ncol(Seurat::Embeddings(obj, "pca")))
  obj <- Seurat::FindNeighbors(obj, dims = 1:nd, verbose = FALSE)
  obj <- Seurat::FindClusters(obj, resolution = 0.4, verbose = FALSE)

  .aa_log(dir, "[prep] clusters computed; K=%d",
          length(unique(obj$seurat_clusters)))
  obj
}

.annot_singleR <- function(obj,
                           species,
                           assay,
                           level = c("cluster", "cell")) {
  level <- match.arg(level)

  feat_ids <- rownames(Seurat::GetAssayData(obj, assay = assay, layer = "counts"))
  use_ensembl <- any(grepl("^ENSG|^ENSMUSG", feat_ids))
  ref <- .pick_celldex(species, use_ensembl = use_ensembl)

  counts <- Seurat::GetAssayData(obj, assay = assay, layer = "counts")

  if (level == "cluster") {
    obj <- .ensure_clusters(obj, assay)
    clusters <- as.character(obj@meta.data$seurat_clusters)
    names(clusters) <- colnames(obj)

    pred <- SingleR::SingleR(
      test = counts,
      ref = ref,
      labels = ref$label.main,
      clusters = clusters,
      assay.type.test = 1
    )

    lab_per_cluster <- setNames(pred$labels, rownames(pred))
    cell_labels <- unname(lab_per_cluster[clusters])
    names(cell_labels) <- colnames(obj)

    conf <- if ("delta.mean" %in% colnames(pred))
      pred$delta.mean[match(clusters, rownames(pred))]
    else
      NULL
    list(cell_labels = cell_labels, conf = conf)

  } else {
    pred <- SingleR::SingleR(
      test = counts,
      ref = ref,
      labels = ref$label.main,
      assay.type.test = 1
    )
    cell_labels <- pred$labels
    names(cell_labels) <- colnames(obj)
    conf <- if ("delta.mean" %in% colnames(pred))
      pred$delta.mean
    else
      NULL
    list(cell_labels = cell_labels, conf = conf)
  }
}

.annot_label_transfer <- function(obj, ref_obj, ref_label_col, assay) {
  stopifnot(ref_label_col %in% colnames(ref_obj@meta.data))
  Seurat::DefaultAssay(ref_obj) <- assay
  ref_obj <- Seurat::NormalizeData(ref_obj, verbose = FALSE)
  ref_obj <- Seurat::FindVariableFeatures(ref_obj, nfeatures = 2000, verbose =
                                            FALSE)
  Seurat::DefaultAssay(obj) <- assay
  obj <- Seurat::NormalizeData(obj, verbose = FALSE)
  obj <- Seurat::FindVariableFeatures(obj, nfeatures = 2000, verbose = FALSE)
  anchors <- Seurat::FindTransferAnchors(reference = ref_obj,
                                         query = obj,
                                         dims = 1:30)
  preds <- Seurat::TransferData(
    anchorset = anchors,
    refdata = ref_obj@meta.data[[ref_label_col]],
    dims = 1:30
  )
  p <- as.character(preds$predicted.id)
  names(p) <- colnames(obj)
  p
}


.ensure_layers_and_data <- function(obj, assay) {
  dir <- .aa_dbg_dir("core")
  Seurat::DefaultAssay(obj) <- assay
  lyr <- try(Seurat::GetAssayData(obj, assay = assay, layer = "data"), silent = FALSE)
  if (inherits(lyr, "try-error") || is.null(lyr) || nrow(lyr) == 0L || ncol(lyr) == 0L) {
    .aa_log(dir, "[prep] NormalizeData(): creating 'data' layer for assay=%s", assay)
    obj <- Seurat::NormalizeData(obj, verbose = FALSE)
  } else {
    .aa_log(dir, "[prep] 'data' layer present (genes=%d cells=%d)", nrow(lyr), ncol(lyr))
  }
  obj
}



.annotate_marker_score <- function(obj,
                                   assay = Seurat::DefaultAssay(obj),
                                   species = c("human", "mouse"),
                                   marker_source = c("panglao", "internal", "custom"),
                                   panglao_file = system.file("extdata",
                                                              "PanglaoDB_markers_27_Mar_2020.tsv",
                                                              package = utils::packageName()),
                                   canonical_only = TRUE,
                                   ui_max = 0.20,
                                   tissue = NULL,
                                   min_genes = 5,
                                   group.by = NULL,
                                   verbose = TRUE) {
  species <- match.arg(species)
  marker_source <- match.arg(marker_source)

  if (marker_source == "panglao") {
    sig <- panglao_signatures(
      tsv = panglao_file,
      species = species,
      canonical_only = canonical_only,
      ui_max = ui_max,
      tissue = tissue,
      min_genes = min_genes
    )
  } else if (marker_source == "internal") {
    sig <- generic_signatures(species = species)
  } else {
    # custom
    stop("Provide your own signatures list for marker_source='custom'")
  }

  obj <- .ensure_layers_and_data(obj, assay = assay)

  genes_obj <- rownames(Seurat::GetAssayData(obj, assay = assay, layer = "data"))
  if (!length(genes_obj)) {
    genes_obj <- rownames(Seurat::GetAssayData(obj, assay = assay, layer = "counts"))
  }
  genes_upper <- toupper(genes_obj)
  map_to_object_case <- function(g) {
    idx <- match(toupper(g), genes_upper)
    idx <- idx[!is.na(idx)]
    if (!length(idx))
      character(0)
    else
      genes_obj[idx]
  }
  sig <- lapply(sig, map_to_object_case)
  sig <- sig[vapply(sig, length, integer(1)) > 0]

  if (!length(sig)) {
    warning("No markers matched object features; skipping marker_score.")
    return(obj)
  }

  # Module scores
  obj <- Seurat::AddModuleScore(
    obj,
    features = sig,
    assay = assay,
    name = "PANG_",
    nbin = 24,
    ctrl = 50
  )
  new_cols <- tail(colnames(obj@meta.data), length(sig))
  new_named <- paste0("PANG_", names(sig))
  for (i in seq_along(new_cols)) {
    j <- match(new_cols[i], colnames(obj@meta.data))
    if (!is.na(j))
      colnames(obj@meta.data)[j] <- new_named[i]
  }
  scmat <- as.matrix(obj@meta.data[, new_named, drop = FALSE])

  obj$auto_celltype      <- sub("^PANG_", "", names(sig))[max.col(scmat, ties.method = "first")]
  obj$auto_celltype_src  <- "marker_score"

  cl <- if (!is.null(group.by) &&
            group.by %in% colnames(obj@meta.data))
    obj[[group.by]][, 1]
  else
    Seurat::Idents(obj)
  maj <- tapply(obj$auto_celltype, cl, function(x)
    names(sort(table(x), TRUE))[1])
  obj$auto_celltype_cluster <- unname(maj[as.character(cl)])
  obj
}


.annot_marker_score_generic <- function(obj,
                                        species = c("mouse", "human"),
                                        assay   = "RNA",
                                        level   = c("cluster", "cell"),
                                        custom_signatures = NULL,
                                        unknown_min_genes  = 2,
                                        unknown_min_margin = 0.15,
                                        unknown_min_top    = 0.10,
                                        nbin = 24,
                                        ctrl = 50,
                                        debug = getOption("scQCenrich.debug", FALSE),
                                        autotune_thresholds = TRUE,
                                        autotune_quantiles = c(top = 0.90, margin = 0.90),
                                        # NEW:
                                        method = c("module_score", "findmarkers"),
                                        panglao_file = system.file("extdata",
                                                                   "PanglaoDB_markers_27_Mar_2020.tsv",
                                                                   package = utils::packageName()),
                                        panglao_tissue = NULL,
                                        panglao_ui_max = 0.20,
                                        panglao_canonical_only = TRUE,
                                        fm_top_n = 150,
                                        fm_padj_max = 1.1,
                                        overlap_mode = c("count", "jaccard")) {
  stopifnot(inherits(obj, "Seurat"))
  species <- match.arg(species)
  level   <- match.arg(level)
  method  <- match.arg(method)
  overlap_mode <- match.arg(overlap_mode)

  dbgdir <- .fm_dbg_dir()
  .log <- function(fmt, ...) .fm_log(dbgdir, fmt, ...)

  .log("[marker_score] begin (method=%s, level=%s)", method, level)

  Seurat::DefaultAssay(obj) <- assay
  lyr <- try(Seurat::GetAssayData(obj, assay = assay, layer = "data"), silent = FALSE)
  if (inherits(lyr, "try-error") || is.null(lyr) || !nrow(lyr) || !ncol(lyr)) {
    .log("[prep] NormalizeData(): creating 'data' layer")
    obj <- Seurat::NormalizeData(obj, verbose = FALSE)
  }
  .log("[prep] computing clusters (assay=%s)", assay)
  obj <- .ensure_clusters(obj, assay = assay)
  K <- length(unique(obj@meta.data$seurat_clusters))
  .log("[prep] clusters computed; K=%d", K)
  .log("[context] cells=%d features=%d", ncol(obj),
       nrow(Seurat::GetAssayData(obj, assay = assay, layer = "data")))

  # ---- PATH A: AddModuleScore (existing behavior) ----
  if (identical(method, "module_score")) {
    # 1) signatures (generic or user)
    sig0 <- if (!is.null(custom_signatures))
      custom_signatures
    else
      generic_signatures(species = species)
    .say("[annot] module_score: %d raw signatures", length(sig0))

    # 2) map genes -> object feature-case (case-insensitive)
    genes_obj   <- rownames(Seurat::GetAssayData(obj, assay = assay, layer = "data"))
    if (!length(genes_obj))
      genes_obj <- rownames(Seurat::GetAssayData(obj, assay = assay, layer = "counts"))
    genes_upper <- toupper(genes_obj)
    map_to_obj  <- function(g) {
      idx <- match(toupper(g), genes_upper)
      idx <- idx[!is.na(idx)]
      if (!length(idx))
        character(0)
      else
        genes_obj[idx]
    }
    sig <- lapply(sig0, map_to_obj)
    keep_mask <- vapply(sig, length, 1L) >= unknown_min_genes
    sig <- sig[keep_mask]
    if (!length(sig))
      return(NULL)

    # 3) module scores
    obj <- Seurat::AddModuleScore(
      obj,
      features = sig,
      assay = assay,
      name = "MS_",
      nbin = nbin,
      ctrl = ctrl
    )
    score_cols <- tail(colnames(obj@meta.data), length(sig))
    colnames(obj@meta.data)[match(score_cols, colnames(obj@meta.data))] <- paste0("MS_", names(sig))
    scmat <- as.matrix(obj@meta.data[, paste0("MS_", names(sig)), drop = FALSE])

    if (level == "cell") {
      top_idx <- max.col(scmat, ties.method = "first")
      top     <- scmat[cbind(seq_len(nrow(scmat)), top_idx)]
      second  <- apply(scmat, 1, function(x)
        .top2(x)[["second"]])
      margin  <- top - second
      labels  <- colnames(scmat)[top_idx]
      labels[(top < unknown_min_top) |
               (margin < unknown_min_margin)] <- "Unknown"
      names(labels) <- rownames(scmat)
      return(list(
        level = "cell",
        labels = labels,
        margin = setNames(margin, rownames(scmat))
      ))
    }

    # cluster-level scoring
    cl  <- as.character(obj$seurat_clusters)
    names(cl) <- colnames(obj)
    df <- data.frame(cluster = cl, scmat, check.names = FALSE)
    cl_means <- stats::aggregate(
      df[, -1, drop = FALSE],
      by = list(cluster = df$cluster),
      FUN = function(x)
        mean(x, na.rm = TRUE)
    )
    rownames(cl_means) <- cl_means$cluster
    cl_means$cluster <- NULL
    cl_means <- as.matrix(cl_means)
    top_idx <- max.col(cl_means, ties.method = "first")
    topC    <- cl_means[cbind(seq_len(nrow(cl_means)), top_idx)]
    secondC <- apply(cl_means, 1, function(x)
      .top2(x)[["second"]])
    marginC <- topC - secondC
    clabels <- colnames(cl_means)[top_idx]
    clabels[(topC < unknown_min_top) |
              (marginC < unknown_min_margin)] <- "Unknown"
    names(clabels) <- rownames(cl_means)

    labels <- clabels[cl]
    names(labels) <- names(cl)
    margin <- setNames(marginC[cl], names(cl))
    return(list(
      level = "cluster",
      labels = labels,
      margin = margin,
      cluster_labels = clabels
    ))
  }

  # ---- PATH B: Panglao overlap of top FindAllMarkers (no AddModuleScore) ----
  .log("[annot] method='findmarkers' using Panglao overlap (top N=%d, padj<=%g, mode=%s)",
       fm_top_n, fm_padj_max, overlap_mode)
  .log("[start] species=%s assay=%s clusters=%d cells=%d features(data)=%d",
       species, assay, K, ncol(obj),
       nrow(Seurat::GetAssayData(obj, assay = assay, layer = "data")))

  # 1) Load Panglao or custom signatures
  sig0 <- if (!is.null(custom_signatures)) {
    .log("custom_signatures provided: %d sets", length(custom_signatures))
    custom_signatures
  } else {
    if (!nzchar(panglao_file) || !file.exists(panglao_file)) {
      .log("[abort] panglao_file missing: %s", panglao_file)
      warning("Panglao TSV not found; provide panglao_file or custom_signatures.")
      return(NULL)
    }
    .log("panglao_signatures(): file=%s species=%s tissue=%s canonical_only=%s ui_max=%g",
         panglao_file, species,
         paste(if (is.null(panglao_tissue)) "<none>" else panglao_tissue, collapse="|"),
         panglao_canonical_only, panglao_ui_max)
    panglao_signatures(
      tsv = panglao_file,
      species = species,
      canonical_only = panglao_canonical_only,
      ui_max = panglao_ui_max,
      tissue = panglao_tissue,
      min_genes = unknown_min_genes
    )
  }
  if (!length(sig0)) {
    .log("[abort] 0 Panglao signatures after filters (tissue/UI/canonical/min_genes)")
    return(NULL)
  }
  .log("panglao_signatures(): kept %d signatures", length(sig0))

  # coverage
  genes_obj <- rownames(Seurat::GetAssayData(obj, assay = assay, layer = "data"))
  if (!length(genes_obj))
    genes_obj <- rownames(Seurat::GetAssayData(obj, assay = assay, layer = "counts"))
  up_obj <- unique(toupper(genes_obj))
  cov_df <- data.frame(
    signature = names(sig0),
    n_markers = vapply(sig0, length, integer(1)),
    overlap   = vapply(sig0, function(v) sum(unique(toupper(v)) %in% up_obj), integer(1)),
    stringsAsFactors = FALSE, check.names = FALSE
  )
  cov_df$frac <- with(cov_df, ifelse(n_markers > 0, overlap / n_markers, 0))
  .fm_csv(dbgdir, "01_signature_coverage.csv", cov_df, row.names = FALSE)
  .log("signature coverage: median overlap = %.2f (IQR %.2f-%.2f)",
       stats::median(cov_df$frac, na.rm = TRUE),
       as.numeric(stats::quantile(cov_df$frac, .25, na.rm = TRUE)),
       as.numeric(stats::quantile(cov_df$frac, .75, na.rm = TRUE)))

  # 2) DE per cluster
  Seurat::Idents(obj) <- "seurat_clusters"
  run_findall <- function(min.pct = 0.10, logfc = 0.25, tag = "strict") {
    de <- try(Seurat::FindAllMarkers(
      obj, assay = assay, only.pos = TRUE,
      logfc.threshold = logfc, min.pct = min.pct, verbose = FALSE
    ), silent = FALSE)
    if (inherits(de, "try-error") || !is.data.frame(de) || !nrow(de)) {
      .log("[FindAllMarkers:%s] failed or 0 rows (min.pct=%.2f logfc=%.2f)", tag, min.pct, logfc)
      return(NULL)
    }
    .log("[FindAllMarkers:%s] rows=%d, cols=%s", tag, nrow(de), paste(colnames(de), collapse=","))
    .fm_csv(dbgdir, sprintf("02_findallmarkers_head_%s.csv", tag), .fm_head_df(de, 25))
    de
  }
  de <- run_findall(0.10, 0.25, "strict")
  if (is.null(de)) de <- run_findall(0.05, 0.10, "relaxed")
  if (is.null(de)) return(NULL)

  # 3) keep columns robustly across Seurat versions
  fc_col   <- intersect(c("avg_log2FC", "avg_logFC"), colnames(de))[1]
  padj_col <- intersect(c("p_val_adj", "p_val_adj_ranked"), colnames(de))[1]
  if (!length(fc_col)) stop("FindAllMarkers output missing logFC column.")
  if (!length(padj_col) && "p_val" %in% colnames(de)) padj_col <- "p_val"  # fallback
  n_before <- nrow(de)
  de <- de[order(de$cluster, -abs(de[[fc_col]]), de[[padj_col]]), , drop=FALSE]
  if (padj_col %in% colnames(de)) {
    keep <- is.na(de[[padj_col]]) | de[[padj_col]] <= fm_padj_max
    de   <- de[keep, , drop = FALSE]
  }
  .log("marker filter: kept %d/%d rows with %s <= %.3g (or NA)", nrow(de), n_before, padj_col, fm_padj_max)

  # 4) select top N per cluster
  de$gene <- as.character(de$gene)
  per_cl  <- stats::aggregate(de$gene, by = list(cluster = de$cluster), FUN = function(v) length(unique(v)))
  names(per_cl)[2] <- "n_unique_genes"
  .fm_csv(dbgdir, "03_markers_per_cluster.csv", per_cl, row.names = FALSE)
  top_by_cl <- split(de$gene, de$cluster)
  top_by_cl <- lapply(top_by_cl, function(v) unique(utils::head(v, fm_top_n)))
  for (k in head(names(top_by_cl), 10)) {
    p <- file.path(dbgdir, sprintf("top_markers_cluster_%s.txt", k))
    writeLines(top_by_cl[[k]], p, useBytes = TRUE)
  }
  .log("top_by_cl: clusters=%d, median(top list size)=%.0f",
       length(top_by_cl), stats::median(vapply(top_by_cl, length, 1L)))

  # 5) overlap scoring
  U <- function(x) unique(toupper(x))
  sigU <- lapply(sig0, U)
  topU <- lapply(top_by_cl, U)
  score_one <- function(top_vec, sig_list, mode = "count") {
    sapply(sig_list, function(ref) {
      if (!length(top_vec) || !length(ref)) return(0)
      inter <- length(intersect(top_vec, ref))
      if (identical(mode, "jaccard")) {
        denom <- length(union(top_vec, ref)); if (denom == 0) 0 else inter / denom
      } else inter
    })
  }
  cl_ids <- names(topU)
  S <- matrix(0, nrow = length(cl_ids), ncol = length(sigU),
              dimnames = list(cl_ids, names(sigU)))
  for (k in cl_ids) S[k, ] <- score_one(topU[[k]], sigU, mode = overlap_mode)
  .fm_csv(dbgdir, "04_overlap_matrix.csv", S)
  ov_summ <- t(apply(S, 1, function(v){o <- sort(v, TRUE); c(top=o[1], second=ifelse(length(o)>=2,o[2], -Inf))}))
  ov_summ <- as.data.frame(ov_summ); ov_summ$cluster <- rownames(ov_summ)
  .fm_csv(dbgdir, "05_overlap_summary_raw.csv", ov_summ)

  # 6) gating + auto-relax if needed
  top_idx <- max.col(S, ties.method = "first")
  top_raw <- S[cbind(seq_len(nrow(S)), top_idx)]
  sec_raw <- apply(S, 1, function(x){o <- sort(x, TRUE); if (length(o)>=2) o[2] else -Inf})
  if (identical(overlap_mode, "count")) {
    top_norm <- top_raw / pmax(fm_top_n, 1); sec_norm <- sec_raw / pmax(fm_top_n, 1)
  } else { top_norm <- top_raw; sec_norm <- sec_raw }
  margin_norm <- top_norm - sec_norm
  clabels <- colnames(S)[top_idx]
  clabels[(top_norm < unknown_min_top) | (margin_norm < unknown_min_margin)] <- "Unknown"
  unk_rate <- mean(clabels == "Unknown")
  .log("gating: min_top=%.3f min_margin=%.3f -> Unknown rate = %.1f%%",
       unknown_min_top, unknown_min_margin, 100 * unk_rate)

  if (all(clabels == "Unknown")) {
    adj <- .auto_relax(top_norm, margin_norm,
                       old_top = unknown_min_top,
                       old_margin = unknown_min_margin,
                       label = "findmarkers")
    clabels2 <- colnames(S)[top_idx]
    clabels2[(top_norm < adj$top) | (margin_norm < adj$margin)] <- "Unknown"
    clabels <- clabels2
    .log("auto-relax applied: new min_top=%.3f min_margin=%.3f -> Unknown rate = %.1f%%",
         adj$top, adj$margin, 100 * mean(clabels == "Unknown"))
  }


  rnS <- rownames(S)
  names(top_norm)    <- rnS
  names(sec_norm)    <- rnS
  names(margin_norm) <- rnS
  names(clabels)     <- rnS

  # --- expand to cells as before ---
  cl  <- as.character(obj$seurat_clusters); names(cl) <- colnames(obj)
  labels_cell <- clabels[cl];           names(labels_cell) <- names(cl)
  margin_cell <- margin_norm[cl];       names(margin_cell) <- names(cl)

  # --- extra debug so you can triage quickly if anything still looks wrong ---
  na_n  <- sum(is.na(labels_cell))
  unk_n <- sum(labels_cell == "Unknown", na.rm = TRUE)
  .log("[postmap] NA labels=%d | Unknown cells=%d (%.1f%%)",
       na_n, unk_n, 100*unk_n/length(labels_cell))
  if (na_n > 0) {
    .log("[warn/postmap] mapping produced NA labels; writing 06_postmap_diag.csv")
    diag_df <- data.frame(
      cell   = names(cl),
      cluster= cl,
      cluster_label   = clabels[cl],
      cluster_margin  = margin_norm[cl],
      stringsAsFactors = FALSE, check.names = FALSE
    )
    .fm_csv(dbgdir, "06_postmap_diag.csv",
            diag_df[seq_len(min(5000, nrow(diag_df))), ], row.names = FALSE)
  }

  # 7) expand to cells; log concise result block
  cl  <- as.character(obj$seurat_clusters); names(cl) <- colnames(obj)
  labels_cell <- clabels[cl]; names(labels_cell) <- names(cl)
  margin_cell <- setNames(margin_norm[cl], names(cl))

  out <- list(level = "cluster",
              labels = labels_cell,
              margin = margin_cell,
              cluster_labels = clabels)
  attr(out, "debug") <- list(
    debug_dir     = dbgdir,
    overlap_mode  = overlap_mode,
    fm_top_n      = fm_top_n,
    padj_col      = padj_col,
    fc_col        = fc_col,
    unknown_gate  = list(min_top = unknown_min_top, min_margin = unknown_min_margin),
    top_norm      = setNames(top_norm, rownames(S)),
    margin_norm   = setNames(margin_norm, rownames(S))
  )
  .log("[done] assigned labels for %d cells across %d clusters",
       length(labels_cell), length(unique(cl)))

  # one-page, copy-friendly result
  snip <- c(
    "[findmarkers] RESULT",
    sprintf(" clusters: %d | cells: %d | overlap_mode: %s | topN: %d | padj_col: %s | fc_col: %s",
            length(unique(cl)), length(labels_cell), overlap_mode, fm_top_n, padj_col, fc_col),
    sprintf(" gate(min_top=%.3f, min_margin=%.3f) -> Unknown=%.1f%%",
            unknown_min_top, unknown_min_margin, 100*mean(clabels == "Unknown")),
    sprintf(" per-cluster(best, top, sec, margin): %s",
            paste(utils::head(sprintf("%s=%s(%.3f,%.3f,%.3f)",
                                      names(clabels), clabels,
                                      as.numeric(top_norm)[match(names(clabels), rownames(S))],
                                      as.numeric(sec_norm)[match(names(clabels), rownames(S))],
                                      as.numeric(margin_norm)[match(names(clabels), rownames(S))]), 6),
                  collapse = " | ")),
    sprintf(" labels(head by cell): %s",
            paste(sprintf("%s=%s", utils::head(names(labels_cell), 6),
                          utils::head(labels_cell, 6)), collapse = ", "))
  )
  # write snippet file and also log it (works even with % in text)
  writeLines(snip, file.path(dbgdir, "RESULT_SNIPPET.txt"))
  .log(paste(snip, collapse = "\n"))

  out
}


# small helper used above to relax gates when everything becomes "Unknown"
.auto_relax <- function(top_vals,
                        margin_vals,
                        old_top,
                        old_margin,
                        label = "cluster") {
  q10_top    <- as.numeric(stats::quantile(top_vals, 0.10, na.rm = TRUE))
  q10_margin <- as.numeric(stats::quantile(margin_vals, 0.10, na.rm = TRUE))
  new_top    <- min(old_top, q10_top)
  new_margin <- min(old_margin, q10_margin)
  message(
    sprintf(
      "[marker_score] AUTO-RELAX(%s): old top=%.3f margin=%.3f -> new top=%.3f margin=%.3f",
      label,
      old_top,
      old_margin,
      new_top,
      new_margin
    )
  )
  list(top = new_top, margin = new_margin)
}

# ---------- ScType backends (sccca multi-sheet or flat Excel) ----------

.annot_sctype <- function(obj,
                          species,
                          assay,
                          level   = c("cluster", "cell"),
                          tissues = NULL,
                          tt      = c("a", "n", "c"),
                          db_path = NULL) {
  level <- match.arg(level)
  tt    <- match.arg(tt)

  if (!requireNamespace("openxlsx", quietly = TRUE)) {
    .dbg("ScType requires 'openxlsx'.")
    return(NULL)
  }

  conv <- try(.to_symbol_assay(obj, species = species, assay = assay),
              silent = F)
  if (!inherits(conv, "try-error")) {
    obj <- conv$obj
    assay <- conv$assay
  }

  obj <- .ensure_clusters(obj, assay = assay)
  lyr_ok <- try(Seurat::GetAssayData(obj, assay = assay, layer = "data"),
                silent = F)
  if (inherits(lyr_ok, "try-error") ||
      nrow(lyr_ok) == 0)
    obj <- Seurat::NormalizeData(obj, verbose = FALSE)

  cls <- as.character(obj$seurat_clusters)
  names(cls) <- colnames(obj)

  tissues <- if (is.null(tissues) ||
                 !length(tissues))
    .sctype_default_tissues()
  else
    as.character(tissues)
  .dbg("[ScType] Using tissues = %s (tt=%s)",
       paste(tissues, collapse = ", "),
       tt)

  if (is.null(db_path) ||
      !nzchar(db_path) || !file.exists(db_path)) {
    .dbg("[WARN] No local ScType Excel provided; skipping ScType backend.")
    return(NULL)
  }

  sh <- try(openxlsx::getSheetNames(db_path), silent = F)
  if (!inherits(sh, "try-error") && length(sh) <= 1L) {
    res <- try(.annot_sctype_flat_excel(
      obj = obj,
      species = species,
      assay = assay,
      db_path = db_path,
      level = level,
      tissues = tissues
    ),
    silent = F)
    if (!inherits(res, "try-error") && !is.null(res))
      return(res)
    return(NULL)
  }

  if (!requireNamespace("sccca", quietly = TRUE)) {
    .dbg("ScType unavailable: install 'sccca'.")
    return(NULL)
  }

  feats <- .cap_features_for_sctype(
    obj,
    assay = assay,
    max_genes = getOption("scQCenrich.sctype_max_genes", 5000L)
  )
  if (length(feats) < 200)
    return(NULL)
  obj <- subset(obj, features = feats)

  cor_m <- try(sccca::calculate_cor_mat(expression_mat = obj,
                                        clusters = "seurat_clusters",
                                        assay = assay),
               silent = F)
  if (inherits(cor_m, "try-error") ||
      is.null(cor_m) || !NROW(cor_m))
    return(NULL)

  org <- if (species == "human")
    "h"
  else
    "m"

  best <- NULL
  for (ti in tissues) {
    sobj2 <- try(sccca::sccca(
      sobj    = obj,
      assay = assay,
      cluster = "seurat_clusters",
      marker  = db_path,
      tissue  = ti,
      tt = tt,
      org     = org,
      test = "p"
    ),
    silent = F)
    if (inherits(sobj2, "try-error") ||
        !inherits(sobj2, "Seurat"))
      next
    labs <- as.character(sobj2@meta.data$cell_type)
    if (length(labs) != ncol(obj))
      next
    names(labs) <- colnames(obj)
    unk <- mean(labs == "Unknown", na.rm = TRUE)
    if (is.null(best) ||
        unk < best$unknown)
      best <- list(tissue = ti,
                   unknown = unk,
                   labels = labs)
  }
  if (!is.null(best) && best$unknown < 0.99) {
    return(list(
      level = level,
      cell_labels = best$labels,
      conf = NULL,
      cluster_labels = NULL
    ))
  }

  # fallback: process_markers + process_clus
  best_lab <- setNames(rep("Unknown", length(unique(cls))), sort(unique(cls)))
  best_p   <- setNames(rep(Inf, length(unique(cls))), sort(unique(cls)))
  for (ti in tissues) {
    db_markers <- try(sccca::process_markers(
      database_path = db_path,
      org = org,
      tissue = ti,
      tissue_type = tt
    ),
    silent = F)
    if (inherits(db_markers, "try-error") ||
        is.null(db_markers))
      next

    for (k in sort(unique(cls))) {
      res <- try(sccca::process_clus(
        cluster = k,
        sobj = obj,
        assay = assay,
        clus = cls,
        markers = db_markers,
        cor_m = cor_m,
        m_t = 0.9,
        c_t = 0.7,
        test = "p"
      ),
      silent = F)
      if (inherits(res, "try-error") ||
          is.null(res) || !NROW(res))
        next
      res <- res[order(res$p_value, -res$overlap_percent, na.last = NA), , drop = FALSE]
      p <- suppressWarnings(min(res$p_value, na.rm = TRUE))
      lab <- res$cell_type[1]
      if (is.finite(p) && p < best_p[[as.character(k)]]) {
        best_p[[as.character(k)]]   <- p
        best_lab[[as.character(k)]] <- paste0(lab, " (", ti, ")")
      }
    }
  }

  if (all(best_lab == "Unknown"))
    return(NULL)

  cell_labels <- unname(best_lab[cls])
  names(cell_labels) <- names(cls)
  conf <- -log10(pmax(best_p[cls], .Machine$double.eps))
  names(conf) <- names(cls)
  list(
    level = "cluster",
    cell_labels = cell_labels,
    conf = conf,
    cluster_labels = best_lab
  )
}

.sctype_default_tissues <- function() {
  c(
    "Immune system",
    "Liver",
    "Lung",
    "Kidney",
    "Brain",
    "Heart",
    "Intestine",
    "Stomach",
    "Pancreas",
    "Skin",
    "Spleen",
    "Thymus",
    "Placenta",
    "Muscle"
  )
}

# flat Excel reader + scorer (tissue-aware)
.read_sctype_flat_xlsx <- function(db_path,
                                   species = c("mouse", "human"),
                                   min_genes = 2,
                                   restrict_tissues = NULL) {
  species <- match.arg(species)
  stopifnot(file.exists(db_path))
  if (!requireNamespace("openxlsx", quietly = TRUE))
    stop("Need openxlsx to read: ", db_path)

  df <- openxlsx::read.xlsx(db_path, sheet = 1)
  if (!is.data.frame(df) ||
      !nrow(df))
    stop("Sheet1 is empty in: ", db_path)

  cn_raw <- names(df)
  cn <- tolower(gsub("\\s+", "_", cn_raw))
  names(df) <- cn

  pick1 <- function(cands) {
    x <- intersect(tolower(cands), names(df))
    if (length(x))
      x[1]
    else
      NA_character_
  }
  col_gene_pos <- pick1(c("genesymbolmore1", "gene_symbol_more_1", "gene"))
  col_gene_neg <- pick1(c("genesymbolmore2", "gene_symbol_more_2"))
  col_ct       <- pick1(c("cellname", "shortname", "cell_type", "celltype", "cell"))
  col_tis      <- pick1(c("tissuetype", "tissue", "organ", "tissue_name"))
  col_sp       <- pick1(c("species", "org", "organism"))

  if (is.na(col_ct) ||
      (is.na(col_gene_pos) && is.na(col_gene_neg))) {
    stop(
      "Need at least a cell-type column and one marker column (geneSymbolmore1/2, cellName)."
    )
  }

  if (!is.na(col_sp)) {
    spv <- tolower(as.character(df[[col_sp]]))
    is_mouse <- grepl("^m(us|m|mouse)|^mm$|^mmu$|^mus musculus$", spv)
    is_human <- grepl("^h(uman|s)?$|^hs$|^hsa$|^homo sapiens$", spv)
    keep <- if (species == "mouse")
      is_mouse |
      (!any(is_mouse, na.rm = TRUE) & !any(is_human, na.rm = TRUE))
    else
      is_human |
      (!any(is_mouse, na.rm = TRUE) & !any(is_human, na.rm = TRUE))
    df <- df[ifelse(is.finite(keep), keep, TRUE), , drop = FALSE]
  }

  if (!is.na(col_tis) && length(restrict_tissues)) {
    tv   <- as.character(df[[col_tis]])
    q    <- unique(na.omit(as.character(restrict_tissues)))
    tv_l <- tolower(trimws(tv))
    q_l  <- tolower(trimws(q))
    keep <- tv_l %in% q_l
    if (!any(keep)) {
      pat  <- paste0(gsub("\\s+", "\\\\s*", q_l), collapse = "|")
      keep <- grepl(pattern = pat,
                    x = tv_l,
                    perl = TRUE)
    }
    if (!any(keep))
      keep <- rep(TRUE, length(tv_l))  # skip filter if nothing matched
    df <- df[keep, , drop = FALSE]
  }

  split_markers <- function(x) {
    if (is.null(x) || length(x) == 0)
      return(character(0))
    x <- as.character(x)
    x[!is.finite(match(x, x))] <- ""
    x <- gsub("///", ",", x, fixed = TRUE)
    parts <- unlist(strsplit(x, ",", fixed = TRUE), use.names = FALSE)
    parts <- toupper(gsub("\\s+", "", parts))
    parts <- parts[parts != "" & parts != "NA"]
    unique(parts)
  }

  ct_vec <- as.character(df[[col_ct]])
  pos <- if (!is.na(col_gene_pos))
    lapply(df[[col_gene_pos]], split_markers)
  else
    replicate(nrow(df), character(0), simplify = FALSE)
  neg <- if (!is.na(col_gene_neg))
    lapply(df[[col_gene_neg]], split_markers)
  else
    replicate(nrow(df), character(0), simplify = FALSE)

  sig_row <- vector("list", length(ct_vec))
  for (i in seq_along(ct_vec)) {
    g <- pos[[i]]
    if (!length(g))
      g <- neg[[i]]
    sig_row[[i]] <- unique(g)
  }

  sig <- tapply(seq_along(ct_vec), ct_vec, function(idx)
    unique(unlist(sig_row[idx], use.names = FALSE)))
  sig <- sig[vapply(sig, length, integer(1)) >= min_genes]

  tissues <- character(0)
  if (!is.na(col_tis) && col_tis %in% names(df)) {
    tv <- as.character(df[[col_tis]])
    tv <- tv[nzchar(tv)]
    if (length(tv))
      tissues <- sort(unique(tv))
  }

  list(signatures = sig, tissues = tissues)
}

.annot_sctype_flat_excel <- function(obj,
                                     species,
                                     assay,
                                     db_path,
                                     level = c("cluster", "cell"),
                                     tissues = NULL,
                                     unknown_min_genes  = 2,
                                     unknown_min_margin = 0.15,
                                     unknown_min_top    = 0.10) {
  level <- match.arg(level)

  parsed <- .read_sctype_flat_xlsx(
    db_path,
    species = species,
    min_genes = unknown_min_genes,
    restrict_tissues = tissues
  )
  sig <- parsed$signatures
  if (!length(sig)) {
    warning("No usable signatures for the requested organ(s).")
    return(NULL)
  }

  .dbg(
    "[ScType] After tissue filter: %d signatures from tissues: %s",
    length(sig),
    paste(unique(tissues), collapse = ", ")
  )

  obj <- .ensure_layers_and_data(obj, assay = assay)

  # map markers -> object feature case
  genes_obj   <- rownames(Seurat::GetAssayData(obj, assay = assay, layer = "data"))
  if (!length(genes_obj))
    genes_obj <- rownames(Seurat::GetAssayData(obj, assay = assay, layer = "counts"))
  genes_upper <- toupper(genes_obj)
  map_to_obj  <- function(g) {
    idx <- match(toupper(g), genes_upper)
    idx <- idx[!is.na(idx)]
    if (!length(idx))
      character(0)
    else
      genes_obj[idx]
  }
  sig <- lapply(sig, map_to_obj)
  sig <- sig[vapply(sig, length, integer(1)) >= unknown_min_genes]
  if (!length(sig)) {
    warning("None of the markers are present in this object after mapping.")
    return(NULL)
  }

  # per-cell module scores
  obj <- Seurat::AddModuleScore(
    obj,
    features = sig,
    assay = assay,
    name = "SCTYPE_",
    nbin = 24,
    ctrl = 50
  )
  new_cols  <- tail(colnames(obj@meta.data), length(sig))
  new_named <- paste0("SCTYPE_", names(sig))
  for (i in seq_along(new_cols)) {
    j <- match(new_cols[i], colnames(obj@meta.data))
    if (!is.na(j))
      colnames(obj@meta.data)[j] <- new_named[i]
  }
  scmat <- as.matrix(obj@meta.data[, new_named, drop = FALSE])

  .top2 <- function(v) {
    o <- sort(v, decreasing = TRUE)
    c(top = o[1], second = ifelse(length(o) >= 2, o[2], -Inf))
  }

  if (level == "cell") {
    top_idx <- max.col(scmat, ties.method = "first")
    top     <- scmat[cbind(seq_len(nrow(scmat)), top_idx)]
    second  <- apply(scmat, 1, function(x)
      .top2(x)[["second"]])
    margin  <- top - second
    labels  <- sub("^SCTYPE_", "", colnames(scmat)[top_idx])
    labels[(top < unknown_min_top) |
             (margin < unknown_min_margin)] <- "Unknown"
    names(labels) <- rownames(scmat)
    return(list(
      level = "cell",
      cell_labels = labels,
      conf = setNames(margin, rownames(scmat))
    ))
  }

  # cluster level
  obj <- .ensure_clusters(obj, assay = assay)
  cl  <- as.character(obj$seurat_clusters)
  names(cl) <- colnames(obj)
  df <- data.frame(cluster = cl, scmat, check.names = FALSE)
  cl_means <- stats::aggregate(
    df[, -1, drop = FALSE],
    by = list(cluster = df$cluster),
    FUN = function(x)
      mean(x, na.rm = TRUE)
  )
  rownames(cl_means) <- cl_means$cluster
  cl_means$cluster <- NULL
  cl_means <- as.matrix(cl_means)

  top_idx <- max.col(cl_means, ties.method = "first")
  top     <- cl_means[cbind(seq_len(nrow(cl_means)), top_idx)]
  second  <- apply(cl_means, 1, function(x)
    .top2(x)[["second"]])
  marginC <- top - second
  clabels <- sub("^SCTYPE_", "", colnames(cl_means)[top_idx])
  names(clabels) <- rownames(cl_means)

  labels <- clabels[cl]
  names(labels) <- names(cl)
  list(
    level = "cluster",
    cell_labels = labels,
    conf = setNames(marginC[cl], names(cl))
  )
}


.write_ms_debug <- function(obj,
                            assay,
                            scmat,
                            sig,
                            top,
                            margin,
                            level,
                            cl_means = NULL,
                            topC = NULL,
                            marginC = NULL,
                            outdir = file.path("qc_outputs", "ms_debug")) {
  dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

  # signature coverage vs object features
  genes_obj <- rownames(Seurat::GetAssayData(obj, assay = assay, layer = "data"))
  if (!length(genes_obj))
    genes_obj <- rownames(Seurat::GetAssayData(obj, assay = assay, layer = "counts"))
  up <- toupper(genes_obj)
  cov <- data.frame(
    signature = names(sig),
    n_markers = vapply(sig, length, integer(1)),
    overlap   = vapply(sig, function(v)
      sum(toupper(v) %in% up), integer(1)),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  cov$frac <- with(cov, ifelse(n_markers > 0, overlap / n_markers, 0))
  utils::write.csv(cov[order(-cov$frac, -cov$overlap), ],
                   file.path(outdir, "01_signature_coverage.csv"),
                   row.names = FALSE)

  # per-cell score head
  utils::write.csv(
    utils::head(as.data.frame(scmat), 5000),
    file.path(outdir, "02_cell_scores_head.csv"),
    row.names = TRUE
  )

  # distributions
  qq <- function(x)
    as.numeric(stats::quantile(x, c(.1, .5, .9), na.rm = TRUE))
  dist <- data.frame(
    metric = c("top", "margin"),
    q10 = c(qq(top)[1], qq(margin)[1]),
    q50 = c(qq(top)[2], qq(margin)[2]),
    q90 = c(qq(top)[3], qq(margin)[3])
  )
  utils::write.csv(dist, file.path(outdir, "03_distributions.csv"), row.names = FALSE)

  if (level == "cluster" && !is.null(cl_means)) {
    utils::write.csv(cl_means,
                     file.path(outdir, "04_cluster_means.csv"),
                     row.names = TRUE)
    cm <- data.frame(cluster = rownames(cl_means),
                     top = topC,
                     margin = marginC)
    utils::write.csv(cm,
                     file.path(outdir, "05_cluster_top_margin.csv"),
                     row.names = FALSE)
  }

  writeLines(paste(names(sig), collapse = ","),
             con = file.path(outdir, "06_signatures_used.txt"))
}

# ---------- Ensembl -> SYMBOL mapper used by ScType ----------

.pick_celldex <- function(species, use_ensembl = FALSE) {
  if (species == "mouse") {
    celldex::MouseRNAseqData(ensembl = use_ensembl)
  } else {
    celldex::HumanPrimaryCellAtlasData(ensembl = use_ensembl)
  }
}

.to_symbol_assay <- function(obj,
                             species = c("mouse", "human"),
                             assay = "RNA",
                             new_assay = "RNA_SYM") {
  species <- match.arg(species)
  rn <- rownames(Seurat::GetAssayData(obj, assay = assay, layer = "counts"))
  looks_ens <- any(grepl("^ENSG\\d+|^ENSMUSG\\d+", rn, ignore.case = FALSE))
  if (!looks_ens)
    return(list(
      obj = obj,
      assay = assay,
      mapped = FALSE
    ))

  if (!requireNamespace("AnnotationDbi", quietly = TRUE))
    stop("Need AnnotationDbi to map Ensembl -> SYMBOL.")
  pkg <- if (species == "mouse")
    "org.Mm.eg.db"
  else
    "org.Hs.eg.db"
  if (!requireNamespace(pkg, quietly = TRUE))
    stop("Need ", pkg, " to map Ensembl -> SYMBOL.")

  suppressWarnings({
    map <- AnnotationDbi::select(
      get(pkg, envir = asNamespace(pkg)),
      keys = rn,
      keytype = "ENSEMBL",
      columns = c("SYMBOL")
    )
  })
  map <- map[!is.na(map$SYMBOL) &
               nzchar(map$SYMBOL), , drop = FALSE]
  if (!nrow(map))
    stop("No Ensembl IDs could be mapped to symbols.")
  map <- map[!duplicated(map$ENSEMBL), , drop = FALSE]
  sym <- setNames(map$SYMBOL, map$ENSEMBL)
  keep <- intersect(names(sym), rn)
  if (!length(keep))
    stop("No overlapping Ensembl IDs after mapping.")

  M <- Seurat::GetAssayData(obj, assay = assay, layer = "counts")
  M <- M[keep, , drop = FALSE]
  grp <- factor(sym[keep], levels = unique(sym[keep]))
  G <- Matrix::sparseMatrix(
    i = as.integer(grp),
    j = seq_along(grp),
    x = 1,
    dims = c(length(levels(grp)), length(grp))
  )
  new_counts <- G %*% M
  rownames(new_counts) <- levels(grp)
  colnames(new_counts) <- colnames(M)

  obj[[new_assay]] <- SeuratObject::CreateAssayObject(counts = new_counts)
  Seurat::DefaultAssay(obj) <- new_assay
  obj <- Seurat::NormalizeData(obj, verbose = FALSE)

  list(obj = obj,
       assay = new_assay,
       mapped = TRUE)
}

.cap_features_for_sctype <- function(obj,
                                     assay = "RNA",
                                     max_genes = 5000) {
  Seurat::DefaultAssay(obj) <- assay
  if (is.null(Seurat::VariableFeatures(obj))) {
    obj <- Seurat::NormalizeData(obj, verbose = FALSE)
    obj <- Seurat::FindVariableFeatures(obj, nfeatures = 8000, verbose = FALSE)
  }
  hvg <- head(Seurat::VariableFeatures(obj), max_genes * 2)
  M   <- Seurat::GetAssayData(obj, assay = assay, layer = "data")[hvg, , drop = FALSE]
  cl  <- as.character(obj$seurat_clusters)
  pb  <- sapply(sort(unique(cl)), function(k)
    Matrix::rowMeans(M[, cl == k, drop = FALSE]))
  keep <- apply(pb, 1, function(v)
    sd(v, na.rm = TRUE)) > 0
  feats <- rownames(M)[keep]
  if (length(feats) > max_genes)
    feats <- feats[seq_len(max_genes)]
  feats
}


.annot_marker_score_debug <- function(obj,
                                      species = c("mouse", "human"),
                                      assay   = "RNA",
                                      level   = c("cluster", "cell"),
                                      signatures = NULL,
                                      unknown_min_genes  = 2,
                                      unknown_min_margin = 0.15,
                                      unknown_min_top    = 0.10,
                                      nbin = 24,
                                      ctrl = 100,
                                      autotune_thresholds = TRUE,
                                      head_n = 6) {
  species <- match.arg(species)
  level   <- match.arg(level)
  # delegate to generic and print a few heads
  res <- .annot_marker_score_generic(
    obj,
    species = species,
    assay = assay,
    level = level,
    custom_signatures = signatures,
    unknown_min_genes = unknown_min_genes,
    unknown_min_margin = unknown_min_margin,
    unknown_min_top = unknown_min_top,
    nbin = nbin,
    ctrl = ctrl,
    debug = TRUE,
    autotune_thresholds = autotune_thresholds
  )
  dbg <- attr(res, "debug")
  if (!is.null(dbg)) {
    message(
      "[MS-DEBUG] level=",
      dbg$level,
      " | n_sig=",
      length(dbg$signature_names),
      " | thresholds: min_top=",
      dbg$thresholds$min_top,
      " min_margin=",
      dbg$thresholds$min_margin
    )
    if (!is.null(dbg$score_head)) {
      print(utils::head(dbg$score_head, head_n))
    }
  }
  res
}

if (!exists(".dbg", mode = "function")) {
  .dbg <- function(fmt, ...) {
    if (isTRUE(getOption("scQCenrich.debug", FALSE))) {
      message(sprintf(fmt, ...))
    }
  }
}



# ---- tiny logging helpers for findmarkers path ----
.fm_dbg_dir <- function() {
  wanted <- file.path("qc_outputs", "findmarkers_debug")
  ok <- dir.create(wanted, recursive = TRUE, showWarnings = FALSE)
  if (!ok)
    wanted <- file.path(tempdir(), "findmarkers_debug")
  dir.create(wanted, recursive = TRUE, showWarnings = FALSE)
  wanted
}
.fm_log <- function(dir, fmt, ...) {
  have_args <- length(list(...)) > 0
  line <- if (have_args) do.call(sprintf, c(list(fmt), list(...))) else fmt
  ts   <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  cat(sprintf("[%s] %s\n", ts, line),
      file = file.path(dir, "log.txt"),
      append = TRUE)
  message(line)  # echo when debug=TRUE
}
.fm_csv <- function(dir, name, x, row.names = TRUE) {
  p <- file.path(dir, name)
  try(utils::write.csv(x, p, row.names = row.names), silent = TRUE)
  return(p)
}
.fm_head_df <- function(x, n = 8) {
  hx <- try(utils::head(as.data.frame(x), n), silent = TRUE)
  if (inherits(hx, "try-error"))
    return(NULL)
  hx
}
.fm_block <- function(title, lines) {
  c(
    paste0("==== ", title, " ===="),
    lines,
    "======================"
  )
}

.annot_log_snippet <- function(obj, res) {
  # build a compact summary without using sprintf
  dbg  <- attr(res, "debug")
  nm   <- names(res$labels)
  labs <- as.character(res$labels)
  unk_rate <- round(100 * mean(labs == "Unknown", na.rm = TRUE), 1)

  # per-cluster summary if available
  cl_summ <- character(0)
  if (!is.null(res$cluster_labels)) {
    tab <- sort(table(res$cluster_labels), decreasing = TRUE)
    cl_summ <- paste0(names(tab), ": ", as.integer(tab))
  }

  head_pairs <- utils::head(paste0(nm, " -> ", labs), 8)

  lines <- c(
    paste0("level: ", res$level),
    paste0("method: ", if (!is.null(dbg$overlap_mode)) paste0("findmarkers/", dbg$overlap_mode) else "module_score"),
    paste0("cells: ", length(labs), " | clusters: ", if (!is.null(res$cluster_labels)) length(unique(names(res$cluster_labels))) else "NA"),
    paste0("Unknown rate: ", unk_rate, "%"),
    if (length(cl_summ)) c("cluster labels (counts):", paste0("  - ", cl_summ)) else NULL,
    "labels (first few):",
    paste0("  - ", head_pairs)
  )
  .fm_block("Auto-annotation summary", lines)
}

# ---- generic debug/log helpers (shared across all backends) ----
`%||%` <- function(x, y) if (is.null(x)) y else x  # small local helper

.aa_dbg_dir <- function(backend = "core") {
  root <- file.path("qc_outputs", "auto_annotate", backend)
  ok <- dir.create(root, recursive = TRUE, showWarnings = FALSE)
  if (!ok) root <- file.path(tempdir(), "auto_annotate", backend)
  dir.create(root, recursive = TRUE, showWarnings = FALSE)
  root
}

.aa_log <- function(dir, fmt, ...) {
  # safe even if there are zero % placeholders
  have_args <- length(list(...)) > 0
  line <- if (have_args) do.call(sprintf, c(list(fmt), list(...))) else fmt
  ts   <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  cat(sprintf("[%s] %s\n", ts, line),
      file = file.path(dir, "log.txt"),
      append = TRUE)
  message(line)
}

.aa_csv <- function(dir, name, x, row.names = TRUE) {
  p <- file.path(dir, name)
  try(utils::write.csv(x, p, row.names = row.names), silent = TRUE)
  p
}

.aa_head_df <- function(x, n = 8) {
  hx <- try(utils::head(as.data.frame(x), n), silent = TRUE)
  if (inherits(hx, "try-error")) return(NULL)
  hx
}

.aa_block <- function(title, lines) {
  c(paste0("==== ", title, " ===="), lines, "======================")
}

.annot_log_snippet <- function(obj, res, backend = "unknown") {
  nm        <- names(res$labels)
  labs      <- as.character(res$labels)
  unk_rate  <- round(100 * mean(labs == "Unknown", na.rm = TRUE), 1)
  cl_summ   <- character(0)

  n_clusters <- "NA"
  if (!is.null(res$cluster_labels)) {
    nn <- names(res$cluster_labels)
    n_clusters <- if (length(nn)) length(unique(nn)) else length(res$cluster_labels)
    tab <- sort(table(res$cluster_labels), decreasing = TRUE)
    cl_summ <- paste0(names(tab), ": ", as.integer(tab))
  }

  head_pairs <- utils::head(paste0(nm, " -> ", labs), 8)
  lines <- c(
    paste0("backend: ", backend),
    paste0("level: ", res$level %||% "<NA>"),
    paste0("cells: ", length(labs), " | clusters: ", n_clusters),
    paste0("Unknown rate: ", unk_rate, "%"),
    if (length(cl_summ)) c("cluster labels (counts):", paste0("  - ", cl_summ)) else NULL,
    "labels (first few):",
    paste0("  - ", head_pairs)
  )
  .aa_block("Auto-annotation summary", lines)
}


