.scq_strip_ver <- function(x) sub("\\.\\d+$", "", x)
.scq_norm_keys <- function(keys, keytype) {
  k <- as.character(keys)
  if (keytype %in% c("ENSEMBL","ENSEMBLTRANS","ENSEMBLPROT","REFSEQ","UNIPROT")) k <- .scq_strip_ver(k)
  if (keytype %in% c("SYMBOL","ALIAS","GENENAME")) k <- toupper(k)
  unique(na.omit(k))
}

.scq_norm_keys_strategy <- function(keys, keytype, strategy=c("norm","raw","title"), organism=c("mouse","human")){
  strategy <- match.arg(strategy); organism <- match.arg(organism)
  k <- as.character(keys)
  .mouse_case <- function(x) sub("^([A-Za-z])([A-Za-z0-9_.-]*)$", "\\U\\1\\L\\2", x, perl = TRUE)
  if (keytype %in% c("ENSEMBL","ENSEMBLTRANS","ENSEMBLPROT","REFSEQ","UNIPROT")) {
    k <- .scq_strip_ver(k)
  } else if (keytype %in% c("SYMBOL","ALIAS","GENENAME")) {
    if (strategy == "norm") {
      k <- toupper(k)
    } else if (strategy == "title") {
      k <- if (organism == "mouse") .mouse_case(k) else k
    } else { # raw
      k <- k
    }
  }
  unique(na.omit(k))
}

.scq_try_map_strategy <- function(keys, OrgDb, keytype,
                                  strategy = c("auto","raw","upper","title")) {
  strategy <- match.arg(strategy)
  # Early exit if OrgDb doesn't support the keytype
  kt_supported <- try(AnnotationDbi::keytypes(OrgDb), silent = TRUE)
  if (inherits(kt_supported, "try-error") || !(keytype %in% kt_supported))
    return(data.frame())

  k <- as.character(keys)
  # strip version for ID-like keytypes
  if (keytype %in% c("ENSEMBL","ENSEMBLTRANS","ENSEMBLPROT","REFSEQ","UNIPROT"))
    k <- sub("\\.\\d+$", "", k)

  # apply case strategy for gene-name-like keytypes
  if (keytype %in% c("SYMBOL","ALIAS","GENENAME")) {
    if (strategy == "upper") {
      k <- toupper(k)
    } else if (strategy == "title") {
      # Mouse-like TitleCase: first letter upper, rest lower (keeps digits/_/.-)
      k <- sub("^([A-Za-z])([A-Za-z0-9_.-]*)$", "\\U\\1\\L\\2", k, perl = TRUE)
    } else if (strategy == "auto") {
      # legacy behavior (kept for regression safety): UPPER
      k <- toupper(k)
    } # "raw" = leave as-is
  }

  k <- unique(stats::na.omit(k))
  if (!length(k)) return(data.frame())

  # Try select(); fall back to mapIds(); never error
  out <- try(AnnotationDbi::select(OrgDb, keys = k, keytype = keytype,
                                   columns = c("ENTREZID","SYMBOL")), silent = TRUE)
  if (inherits(out, "try-error") || !is.data.frame(out)) {
    out <- try({
      ent <- AnnotationDbi::mapIds(OrgDb, keys = k, keytype = keytype,
                                   column = "ENTREZID", multiVals = "first")
      sy  <- suppressWarnings(AnnotationDbi::mapIds(OrgDb, keys = names(ent),
                                                    keytype = keytype,
                                                    column = "SYMBOL",
                                                    multiVals = "first"))
      data.frame(ENTREZID = as.character(ent),
                 SYMBOL   = as.character(sy),
                 key      = names(ent),
                 stringsAsFactors = FALSE)
    }, silent = TRUE)
    if (inherits(out, "try-error") || !is.data.frame(out)) return(data.frame())
  }

  if (!("ENTREZID" %in% names(out))) return(data.frame())
  out$ENTREZID <- as.character(out$ENTREZID)
  out <- out[!is.na(out$ENTREZID) & nzchar(out$ENTREZID), , drop = FALSE]
  unique(out)
}

.scq_try_map <- function(keys, OrgDb, keytype) {
  kt_supported <- try(AnnotationDbi::keytypes(OrgDb), silent = TRUE)
  if (inherits(kt_supported, "try-error") || !(keytype %in% kt_supported)) return(data.frame())
  k <- .scq_norm_keys(keys, keytype)
  if (!length(k)) return(data.frame())
  out <- try(AnnotationDbi::select(OrgDb, keys = unique(k), keytype = keytype,
                                   columns = c("ENTREZID","SYMBOL")), silent = TRUE)
  if (inherits(out, "try-error") || !is.data.frame(out)) {
    out <- try({
      ent <- AnnotationDbi::mapIds(OrgDb, keys = unique(k), keytype = keytype,
                                   column = "ENTREZID", multiVals = "first")
      sy  <- suppressWarnings(AnnotationDbi::mapIds(OrgDb, keys = names(ent), keytype = keytype,
                                                    column = "SYMBOL",  multiVals = "first"))
      data.frame(ENTREZID = as.character(ent), SYMBOL = as.character(sy),
                 KEY = names(ent), stringsAsFactors = FALSE)
    }, silent = TRUE)
    if (inherits(out, "try-error") || !is.data.frame(out)) return(data.frame())
    if (!("ENTREZID" %in% names(out))) return(data.frame())
  }
  out$ENTREZID <- as.character(out$ENTREZID)
  out <- out[!is.na(out$ENTREZID) & nzchar(out$ENTREZID), , drop = FALSE]
  unique(out)
}

.scq_guess_best_keytype <- function(keys, OrgDb, try_order = NULL, sample_n = 5000L) {
  if (is.null(try_order)) try_order <- c("SYMBOL","ALIAS","ENSEMBL","ENSEMBLTRANS","ENSEMBLPROT","REFSEQ","UNIPROT","GENENAME")
  kt_supported <- try(AnnotationDbi::keytypes(OrgDb), silent = TRUE)
  if (inherits(kt_supported, "try-error")) kt_supported <- character(0)
  try_order <- intersect(try_order, kt_supported)
  if (!length(try_order)) return(list(best = "SYMBOL", scores = c(SYMBOL = 0)))
  kk <- as.character(keys); if (length(kk) > sample_n) { set.seed(1L); kk <- sample(kk, sample_n) }
  scores <- sapply(try_order, function(kt) {
    res <- try(.scq_try_map(kk, OrgDb, kt), silent = TRUE)
    if (inherits(res, "try-error") || !nrow(res)) 0L else nrow(res)
  })
  kt_best <- try_order[which.max(scores)]
  list(best = kt_best, scores = sort(scores, decreasing = TRUE))
}


.scq_ensure_seurat_clusters <- function(
    obj,
    assay      = Seurat::DefaultAssay(obj),
    resolution = 0.4,
    dims       = 1:20,
    npcs       = 30
){
  stopifnot(inherits(obj, "Seurat"))
  if ("seurat_clusters" %in% colnames(obj@meta.data)) {
    message(sprintf("[prep] using existing seurat_clusters (K=%d)", length(unique(obj$seurat_clusters))))
    return(obj)
  }

  message(sprintf("[prep] computing seurat_clusters (resolution=%.2f)", resolution))

  id_backup <- Seurat::Idents(obj)
  on.exit({try(Seurat::Idents(obj) <- id_backup, silent = TRUE)}, add = TRUE)

  Seurat::DefaultAssay(obj) <- assay

  if (!length(Seurat::VariableFeatures(obj))) {
    obj <- Seurat::FindVariableFeatures(obj, nfeatures = 2000, verbose = FALSE)
  }

  need_pcs <- max(max(dims), npcs)
  has_pca  <- "pca" %in% names(obj@reductions)
  have_enough_pcs <- has_pca && ncol(Seurat::Embeddings(obj[["pca"]])) >= need_pcs
  if (!have_enough_pcs) {
    obj <- Seurat::ScaleData(obj, verbose = FALSE)
    obj <- Seurat::RunPCA(obj, npcs = need_pcs, verbose = FALSE)
  }

  obj <- Seurat::FindNeighbors(obj, dims = dims, verbose = FALSE)
  obj <- Seurat::FindClusters(obj, resolution = resolution, verbose = FALSE)

  message(sprintf("[prep] seurat_clusters created (K=%d)", length(unique(obj$seurat_clusters))))
  obj
}

.sqe_pick_orgdb <- function(species = c("human","mouse")) {
  species <- match.arg(tolower(species))
  if (species == "human") {
    if (!requireNamespace("org.Hs.eg.db", quietly = TRUE))
      stop("org.Hs.eg.db required but not installed")
    return(org.Hs.eg.db::org.Hs.eg.db)
  } else {
    if (!requireNamespace("org.Mm.eg.db", quietly = TRUE))
      stop("org.Mm.eg.db required but not installed")
    return(org.Mm.eg.db::org.Mm.eg.db)
  }
}

.sqe_norm_keys <- function(x, keytype) {
  x <- as.character(x)
  if (keytype %in% c("ENSEMBL","ENSEMBLTRANS","ENSEMBLPROT","REFSEQ","UNIPROT"))
    x <- sub("\\.\\d+$", "", x)   # strip version suffix (e.g., ENSG... .12)
  if (keytype %in% c("SYMBOL","ALIAS","GENENAME"))
    x <- toupper(x)               # robust case-insensitive match
  unique(stats::na.omit(x))
}

.sqe_bitr_safe <- function(genes, fromType = "SYMBOL",
                           toType = "ENTREZID", OrgDb) {
  genes <- .sqe_norm_keys(genes, fromType)
  if (!length(genes)) {
    return(data.frame(FROM_ID=character(), TO_ID=character())[FALSE,])
  }
  out <- tryCatch(
    suppressWarnings(
      clusterProfiler::bitr(genes, fromType = fromType,
                            toType   = toType,
                            OrgDb    = OrgDb,
                            drop     = TRUE)  # normal to drop some IDs
    ),
    error = function(e) {
      data.frame(FROM_ID=character(), TO_ID=character())[FALSE,]
    }
  )
  if (nrow(out)) {
    names(out)[names(out) == fromType] <- "FROM_ID"
    names(out)[names(out) == toType]   <- "TO_ID"
    out <- unique(out[, c("FROM_ID","TO_ID")])
  }
  out
}
# -----------------------------------------------------------------------------


validation_plots_post_annotation <- function(
    obj,
    group.by = "auto_celltype",
    organism = c("mouse","human"),
    assay = Seurat::DefaultAssay(obj),
    slot.use = "data",
    imgDir = "qc_outputs",
    top_n = 100,
    min_genes_for_enrich = 10,
    top_go_per_cluster = 1
){
  organism <- match.arg(organism)
  if (!dir.exists(imgDir)) dir.create(imgDir, recursive = TRUE, showWarnings = FALSE)
  stopifnot(inherits(obj, "Seurat"))
  if (!(group.by %in% colnames(obj@meta.data))) {
    warning("validation_plots_post_annotation(): '", group.by, "' not found; skipping validation plots.")
    return(invisible(NULL))
  }

  # ---- packages we (may) use (unchanged)
  has_CH <- requireNamespace("ComplexHeatmap", quietly = TRUE) &&
    requireNamespace("circlize", quietly = TRUE) &&
    requireNamespace("grid", quietly = TRUE)
  has_CP <- requireNamespace("clusterProfiler", quietly = TRUE) &&
    requireNamespace("AnnotationDbi", quietly = TRUE)

  OrgDbPkg <- if (organism == "mouse") "org.Mm.eg.db" else "org.Hs.eg.db"
  has_Org  <- requireNamespace(OrgDbPkg, quietly = TRUE)
  OrgDb    <- if (has_Org) getExportedValue(OrgDbPkg, OrgDbPkg) else NULL

  # ---- 1) Markers per group (top features) (unchanged)
  Seurat::Idents(obj) <- obj@meta.data[[group.by]]
  markers_all <- Seurat::FindAllMarkers(
    obj, assay = assay, slot = slot.use,
    only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25, verbose = FALSE
  )
  if (!nrow(markers_all)) {
    warning("No markers found; skipping validation plots.")
    return(invisible(NULL))
  }
  lfcol <- if ("avg_log2FC" %in% colnames(markers_all)) "avg_log2FC" else
    if ("avg_logFC"   %in% colnames(markers_all)) "avg_logFC"   else "p_val_adj"

  # ---- 2) Heatmap gene set (mostly unchanged)
  top_per_cluster10 <- dplyr::group_by(markers_all, cluster) |>
    dplyr::slice_max(order_by = .data[[lfcol]], n = 10, with_ties = FALSE) |>
    dplyr::ungroup()

  genes_use <- unique(top_per_cluster10$gene)

  expr_avg <- try(
    Seurat::AverageExpression(obj, assays = assay, slot = slot.use, features = genes_use, return.seurat = FALSE),
    silent = TRUE
  )
  if (!inherits(expr_avg, "try-error") && is.list(expr_avg) && length(expr_avg) >= 1) {
    mat <- expr_avg[[assay]]
  } else {
    agg <- Seurat::AggregateExpression(obj, assays = assay, slot = slot.use, features = genes_use)
    mat <- agg[[assay]]
  }

  mat <- as.matrix(mat[intersect(rownames(mat), genes_use), , drop = FALSE])
  if (!is.null(rownames(mat))) rownames(mat) <- make.unique(as.character(rownames(mat)))
  norm <- function(x) gsub("-", "_", x, fixed = TRUE)
  colnames(mat) <- norm(colnames(mat))
  tpc_cluster <- norm(top_per_cluster10$cluster)

  md <- obj@meta.data
  all_groups <- colnames(mat)
  cell_counts_tab <- table(factor(norm(md[[group.by]]), levels = all_groups))
  col_order <- names(cell_counts_tab)
  mat <- mat[, col_order[col_order %in% colnames(mat)], drop = FALSE]

  ordered_genes <- character(0)
  for (cl in colnames(mat)) {
    g_this <- top_per_cluster10$gene[tpc_cluster == cl]
    g_this <- intersect(g_this, rownames(mat))
    if (length(g_this)) {
      ord <- order(mat[g_this, cl], decreasing = TRUE, na.last = TRUE)
      ordered_genes <- c(ordered_genes, g_this[ord])
    } else {
      ord <- order(rowMeans(abs(mat[, cl, drop = FALSE]), na.rm = TRUE), decreasing = TRUE)
      take <- head(rownames(mat)[ord], 5)
      ordered_genes <- c(ordered_genes, take)
    }
  }
  ordered_genes <- unique(ordered_genes)
  if (length(ordered_genes) == 0) {
    ord_global <- order(rowMeans(abs(mat), na.rm = TRUE), decreasing = TRUE)
    ordered_genes <- head(rownames(mat)[ord_global], max(1, min(30, nrow(mat))))
  }
  mat <- mat[ordered_genes, , drop = FALSE]

  # ---- 3) Heatmap output (unchanged)
  heat_png <- file.path(imgDir, "celltype_gene_heatmap.png")
  if (has_CH && nrow(mat) && ncol(mat)) {
    finite_vals <- mat[is.finite(mat)]
    qq <- stats::quantile(finite_vals, probs = c(0.02, 0.98), na.rm = TRUE)
    lo <- as.numeric(qq[1]); hi <- as.numeric(qq[2])
    if (!is.finite(lo)) lo <- min(finite_vals, na.rm = TRUE)
    if (!is.finite(hi)) hi <- max(finite_vals, na.rm = TRUE)
    col_fun <- circlize::colorRamp2(seq(lo, hi, length.out = 7),
                                    grDevices::colorRampPalette(RColorBrewer::brewer.pal(7, "YlGnBu"))(7))
    w <- min(2400, 300 + 90 * ncol(mat))
    h <- min(3200, 300 + 22 * nrow(mat))
    grDevices::png(heat_png, width = w, height = h, res = 150)
    ComplexHeatmap::draw(
      ComplexHeatmap::Heatmap(
        mat, name = "expr", col = col_fun,
        show_row_names = TRUE, row_names_gp = grid::gpar(fontsize = 8),
        show_column_names = TRUE, column_names_rot = 90,
        cluster_rows = FALSE, cluster_columns = FALSE,
        row_title = NULL, column_title = group.by
      ),
      heatmap_legend_side = "right"
    )
    grDevices::dev.off()
  } else {
    z <- t(scale(t(mat))); z[!is.finite(z)] <- 0
    zfinite <- z[is.finite(z)]
    zlim <- if (!length(zfinite)) c(0,1) else {
      zr <- range(zfinite, na.rm = TRUE)
      if (!all(is.finite(zr))) c(0,1) else if (zr[1] == zr[2]) {
        eps <- if (zr[1] == 0) 1e-6 else abs(zr[1]) * 1e-6
        c(zr[1]-eps, zr[2]+eps)
      } else zr
    }
    grDevices::png(heat_png, width = 1400, height = 1000, res = 150)
    op <- graphics::par(mar = c(6, 8, 2, 2))
    nx <- ncol(z); ny <- nrow(z)
    graphics::image(x = seq_len(nx), y = seq_len(ny),
                    z = t(z[ny:1, , drop = FALSE]),
                    xlab = "", ylab = "", axes = FALSE, useRaster = TRUE, zlim = zlim)
    graphics::axis(2, at = seq_len(ny), labels = rownames(z)[ny:1], las = 2, cex.axis = 0.6)
    graphics::axis(1, at = seq_len(nx), labels = colnames(z), las = 2, cex.axis = 0.7)
    graphics::par(op); grDevices::dev.off()
  }

  go_png <- file.path(imgDir, "celltype_go_barplot.png")
  go_csv <- file.path(imgDir, "celltype_go_top.csv")
  dbg_go <- file.path(imgDir, "annot_validation", "go")
  dir.create(dbg_go, recursive = TRUE, showWarnings = FALSE)

  top_per_cluster <- dplyr::group_by(markers_all, cluster) |>
    dplyr::slice_max(order_by = .data[[lfcol]], n = top_n, with_ties = FALSE) |>
    dplyr::ungroup()

  keep_nonsig_placeholder   <- TRUE
  export_full_enrich_table  <- TRUE
  export_inputs             <- TRUE
  min_show_unadjusted       <- 1

  if (has_CP && has_Org) {
    feats <- rownames(Seurat::GetAssayData(obj, assay = assay, layer = "counts"))

    # candidate keytypes present in OrgDb
    kt_all <- try(AnnotationDbi::keytypes(OrgDb), silent = TRUE)
    if (inherits(kt_all, "try-error")) kt_all <- character(0)
    kt_cand <- intersect(c("ENTREZID","SYMBOL","ALIAS","ENSEMBL","ENSEMBLTRANS",
                           "ENSEMBLPROT","REFSEQ","UNIPROT","GENENAME"), kt_all)
    if (!length(kt_cand)) kt_cand <- "SYMBOL"

    probe_rows <- list()
    for (kt in kt_cand) {
      for (st in c("raw","title","upper","auto")) {
        res <- try(.scq_try_map_strategy(feats, OrgDb, kt, strategy = st), silent = TRUE)
        nm  <- if (inherits(res, "try-error") || !nrow(res)) 0L else
          length(unique(stats::na.omit(res$ENTREZID)))
        probe_rows[[paste(kt,st,sep=":")]] <- data.frame(KEYTYPE=kt, STRATEGY=st,
                                                         N_MAPPED=nm, stringsAsFactors=FALSE)
      }
    }
    probe_df <- do.call(rbind, probe_rows)
    utils::write.csv(probe_df, file.path(dbg_go, "00_probe_keytype_yields.csv"), row.names = FALSE)

    # pick best (max mappings)
    best_idx <- which.max(probe_df$N_MAPPED)
    id_from  <- probe_df$KEYTYPE[best_idx]
    strategy <- probe_df$STRATEGY[best_idx]
    best_n   <- probe_df$N_MAPPED[best_idx]
    frac     <- ifelse(length(feats), best_n/length(unique(feats)), 0)
    writeLines(sprintf("TOP5:\n%s",
                       paste(utils::head(probe_df[order(-probe_df$N_MAPPED),], 5L)
                             |> apply(1, function(r) paste(r, collapse=" | ")),
                             collapse="\n")),
               file.path(dbg_go, "00_probe_top5.txt"))

    message(sprintf("[GO] Auto-detected keytype: %s | strategy=%s | mapped=%d / %d (%.1f%%)",
                    id_from, strategy, best_n, length(unique(feats)), 100*frac))
    utils::write.csv(probe_df[best_idx, , drop = FALSE],
                     file.path(dbg_go, "00_keytype_guess.csv"), row.names = FALSE)

    # map universe with chosen strategy; if still sparse, hybrid-union rescue
    map_univ <- .scq_try_map_strategy(feats, OrgDb, id_from, strategy = strategy)
    common_cols <- intersect(c("ENTREZID","SYMBOL","ALIAS","ENSEMBL","ENSEMBLTRANS","ENSEMBLPROT"), colnames(map_univ))
    map_univ <- unique(map_univ[, common_cols, drop = FALSE])
    universe_entrez <- unique(stats::na.omit(map_univ$ENTREZID))

    hybrid_note <- NA_character_
    if (length(universe_entrez) < 0.05 * length(unique(feats))) {
      hybrid_note <- "hybrid_union_applied"
      union_maps <- list(map_univ)
      for (kt in intersect(c("SYMBOL","ALIAS","ENSEMBL","ENSEMBLTRANS","ENSEMBLPROT"), kt_cand)) {
        for (st in c("raw","title","upper","auto")) {
          tmp <- .scq_try_map_strategy(feats, OrgDb, kt, strategy = st)
          if (is.data.frame(tmp) && nrow(tmp)) {
            tmp <- tmp[, intersect(c("ENTREZID","SYMBOL","ALIAS","ENSEMBL","ENSEMBLTRANS","ENSEMBLPROT"), colnames(tmp)), drop = FALSE]
            union_maps[[paste(kt,st,sep=":")]] <- tmp
          }
        }
      }
      union_df <- do.call(rbind, union_maps)
      union_df <- unique(union_df[!is.na(union_df$ENTREZID) & nzchar(union_df$ENTREZID), , drop = FALSE])
      universe_entrez <- unique(union_df$ENTREZID)
      utils::write.csv(utils::head(union_df, 100),
                       file.path(dbg_go, "01_universe_union_head.csv"), row.names = FALSE)
    }

    writeLines(c(
      sprintf("features=%d", length(unique(feats))),
      sprintf("universe_mapped=%d", length(universe_entrez)),
      sprintf("keytype=%s", id_from),
      sprintf("strategy=%s", strategy),
      sprintf("hybrid=%s", ifelse(is.na(hybrid_note), "no", hybrid_note))
    ), con = file.path(dbg_go, "01_universe_summary.txt"))
    utils::write.csv(utils::head(map_univ, 50), file.path(dbg_go, "01_universe_head.csv"), row.names = FALSE)

    # map top genes pool using the SAME strategy; if weak, union fallback
    genes_pool <- unique(top_per_cluster$gene)
    map_top <- .scq_try_map_strategy(genes_pool, OrgDb, id_from, strategy = strategy)
    if (!nrow(map_top) || sum(!is.na(map_top$ENTREZID)) < 3) {
      union_top <- list(map_top)
      for (kt in intersect(c("SYMBOL","ALIAS","ENSEMBL","ENSEMBLTRANS","ENSEMBLPROT"), kt_cand)) {
        for (st in c("raw","title","upper","auto")) {
          tmp <- .scq_try_map_strategy(genes_pool, OrgDb, kt, strategy = st)
          if (is.data.frame(tmp) && nrow(tmp))
            union_top[[paste(kt,st,sep=":")]] <- tmp
        }
      }
      map_top <- unique(do.call(rbind, lapply(union_top, function(x){
        x[, intersect(c("ENTREZID","SYMBOL","ALIAS","ENSEMBL","ENSEMBLTRANS","ENSEMBLPROT"), colnames(x)), drop = FALSE]
      })))
    }
    utils::write.csv(utils::head(map_top, 50), file.path(dbg_go, "02_map_top_head.csv"), row.names = FALSE)

    # build key->entrez dictionary using the best column present
    key_col <- intersect(c(id_from,"SYMBOL","ALIAS","ENSEMBL","ENSEMBLTRANS","ENSEMBLPROT"), colnames(map_top))
    key_col <- if (length(key_col)) key_col[1] else "SYMBOL"
    key2ent <- setNames(map_top$ENTREZID, as.character(map_top[[key_col]]))

    # dump unmapped examples for inspection
    unm <- setdiff(genes_pool, unique(as.character(map_top[[key_col]])))
    utils::write.csv(data.frame(unmapped_example = utils::head(unm, 60)),
                     file.path(dbg_go, "02_unmapped_examples.csv"), row.names = FALSE)

    # enrichment per cell type (diagnostics intact + expanded)
    go_results <- list()
    clusters <- colnames(mat)
    min_needed <- min_genes_for_enrich
    if (min_needed > 5 && length(universe_entrez) < 2000) min_needed <- max(5, floor(min_needed * 0.6))

    .diag_template <- function() {
      data.frame(
        celltype   = NA_character_,
        n_top      = NA_integer_,
        n_mapped   = NA_integer_,
        min_needed = NA_integer_,
        universe   = length(universe_entrez),
        status     = NA_character_,
        n_terms_all= NA_integer_,
        n_terms_sig= NA_integer_,
        pmin       = NA_real_,
        padj_min   = NA_real_,
        Count_max  = NA_integer_,
        note       = NA_character_,
        stringsAsFactors = FALSE
      )
    }
    diag_rows <- list()

    for (ct in clusters) {
      g_keys <- unique(top_per_cluster$gene[top_per_cluster$cluster == ct])
      g_entrez <- unique(stats::na.omit(key2ent[as.character(g_keys)]))

      if (export_inputs) {
        df_in <- data.frame(input_key = g_keys,
                            matched_key = names(key2ent)[match(as.character(g_keys), names(key2ent))],
                            ENTREZID = unname(key2ent[as.character(g_keys)]),
                            stringsAsFactors = FALSE)
        # utils::write.csv(df_in, file.path(dbg_go, sprintf("in_%s_genes.csv", ct)), row.names = FALSE)
      }

      message(sprintf("[GO] %s: top=%d -> mapped=%d (min=%d)",
                      ct, length(g_keys), length(g_entrez), min_needed))

      row <- .diag_template()
      row$celltype   <- ct
      row$n_top      <- length(g_keys)
      row$n_mapped   <- length(g_entrez)
      row$min_needed <- min_needed

      if (length(g_entrez) < min_needed || length(universe_entrez) < 100) {
        row$status <- "skipped:too_few_ids_or_universe"
        diag_rows[[length(diag_rows)+1]] <- row
        go_results[[ct]] <- if (keep_nonsig_placeholder) data.frame(
          ID = NA, Description = "(no significant BP term)", p.adjust = 1, Count = 0,
          cluster = ct, stringsAsFactors = FALSE) else NULL
        next
      }

      enr <- tryCatch(
        clusterProfiler::enrichGO(
          gene          = g_entrez,
          OrgDb         = OrgDb,
          keyType       = "ENTREZID",
          ont           = "BP",
          pAdjustMethod = "BH",
          universe      = universe_entrez,
          qvalueCutoff  = 0.2,
          pvalueCutoff  = 0.05,
          readable      = TRUE
        ),
        error = function(e) { message(sprintf("[GO] %s: enrichGO ERROR: %s", ct, e$message)); NULL }
      )

      if (is.null(enr)) {
        row$status <- "skipped:enrich_error"
        diag_rows[[length(diag_rows)+1]] <- row
        go_results[[ct]] <- NULL
        next
      }

      df_full <- as.data.frame(enr)
      n_all <- NROW(df_full); n_sig <- if (n_all) sum(df_full$p.adjust <= 0.20, na.rm = TRUE) else 0L
      # if (export_full_enrich_table) utils::write.csv(df_full, file.path(dbg_go, sprintf("ego_%s_full.csv", ct)), row.names = FALSE)

      if (!n_all) {
        row$status      <- "no_terms"
        row$n_terms_all <- 0L; row$n_terms_sig <- 0L
        diag_rows[[length(diag_rows)+1]] <- row
        go_results[[ct]] <- if (keep_nonsig_placeholder) data.frame(
          ID = NA, Description = "(no significant BP term)", p.adjust = 1, Count = 0,
          cluster = ct, stringsAsFactors = FALSE) else NULL
        next
      }

      row$status      <- if (n_sig) "ok" else "no_sig_terms_after_cutoffs"
      row$n_terms_all <- n_all
      row$n_terms_sig <- n_sig
      row$pmin        <- min(df_full$pvalue,  na.rm=TRUE)
      row$padj_min    <- min(df_full$p.adjust,na.rm=TRUE)
      row$Count_max   <- max(df_full$Count,   na.rm=TRUE)
      diag_rows[[length(diag_rows)+1]] <- row

      if (n_sig == 0 && min_show_unadjusted > 0) {
        df_full <- df_full[order(df_full$pvalue, -df_full$Count), , drop = FALSE]
        df_full <- utils::head(df_full, min_show_unadjusted)
      } else {
        df_full <- df_full[df_full$p.adjust <= 0.20, , drop = FALSE]
      }

      if (NROW(df_full) > 0) {
        df_full$cluster <- rep(ct, NROW(df_full))  # <-- SAFE ADD
        go_results[[ct]] <- df_full
      } else {
        go_results[[ct]] <- if (keep_nonsig_placeholder) data.frame(
          ID = NA, Description = "(no significant BP term)", p.adjust = 1, Count = 0,
          cluster = ct, stringsAsFactors = FALSE) else NULL
      }
    }

    if (length(diag_rows)) {
      diag_df <- do.call(rbind, diag_rows)
      utils::write.csv(diag_df, file.path(dbg_go, "03_per_celltype_diagnostics.csv"), row.names = FALSE)
    }

    # combine, write, plot (unchanged)
    go_df_list <- lapply(names(go_results), function(ct) {
      df <- go_results[[ct]]
      if (is.null(df)) return(NULL)
      df <- df[order(df$p.adjust, -df$Count), ]
      utils::head(df, top_go_per_cluster)[, c("cluster","ID","Description","p.adjust","Count")]
    })
    go_df <- do.call(rbind, go_df_list)

    if (!is.null(go_df) && nrow(go_df)) {
      go_df <- as.data.frame(go_df, stringsAsFactors = FALSE)
      names(go_df)[names(go_df) == "cluster"]     <- "celltype"
      names(go_df)[names(go_df) == "Description"] <- "term"
      go_df$p_neglog10 <- -log10(go_df$p.adjust)
      go_df$celltype <- factor(go_df$celltype, levels = colnames(mat))
      go_df$placeholder <- is.na(go_df$ID) | grepl("^\\(no significant", go_df$term)

      ord <- order(as.integer(go_df$celltype), -go_df$p_neglog10, -go_df$Count)
      go_df <- go_df[ord, , drop = FALSE]
      go_df$term_label <- factor(paste0(go_df$celltype, ": ", go_df$term),
                                 levels = rev(unique(paste0(go_df$celltype, ": ", go_df$term))))
      utils::write.csv(go_df, go_csv, row.names = FALSE)

      if (requireNamespace("ggplot2", quietly = TRUE)) {
        gg <- ggplot2::ggplot(go_df, ggplot2::aes(x = p_neglog10, y = term_label, fill = celltype, alpha = !placeholder)) +
          ggplot2::geom_col() +
          ggplot2::scale_alpha_manual(values = c(`TRUE` = 0.4, `FALSE` = 1), guide = "none") +
          ggplot2::labs(x = expression(-log[10]("adj. P")), y = NULL, fill = "Cell type",
                        title = "Top GO BP term per cell type") +
          ggplot2::theme_minimal(base_size = 10) +
          ggplot2::theme(legend.position = "none",
                         panel.grid = ggplot2::element_blank(),
                         axis.text.y = ggplot2::element_text(size = 9))
        ggplot2::ggsave(go_png, gg, width = 8.5, height = nrow(go_df)*0.15 + 1,
                        dpi = 150, bg = "transparent")
      }
    } else {
      utils::write.csv(data.frame(), go_csv, row.names = FALSE)
      writeLines(c(
        "==== GO enrichment skipped ====",
        sprintf("Universe mapped: %d", length(universe_entrez)),
        sprintf("Chosen keytype: %s (strategy=%s)", id_from, strategy),
        sprintf("Hybrid: %s", ifelse(is.na(hybrid_note), "no", hybrid_note)),
        "See 00_probe_keytype_yields.csv, 00_probe_top5.txt and 03_per_celltype_diagnostics.csv"
      ), file.path(dbg_go, "FAIL_SNIPPET.txt"))
    }
  } else {
    warning("clusterProfiler/AnnotationDbi/OrgDb not available; GO barplot skipped. Install org.*.eg.db for your organism.")
  }

  # --- 4b) GO BP per Seurat cluster (robust + verbose; always writes a PNG)
  cluster_go_png <- file.path(imgDir, "cluster_go_barplot.png")
  cluster_go_csv <- file.path(imgDir, "cluster_go_top_per_cluster.csv")
  obj <- .scq_ensure_seurat_clusters(obj, assay = assay, resolution = 0.4, dims = 1:20, npcs = 30)

  if (requireNamespace("clusterProfiler", quietly = TRUE) &&
      requireNamespace("AnnotationDbi", quietly = TRUE) &&
      !is.null(OrgDb) && inherits(OrgDb, "OrgDb")) {

    # Majority cell-type label per Seurat cluster (for y-axis label)
    cl_vec <- as.character(obj$seurat_clusters)
    ct_vec <- as.character(obj@meta.data[[group.by]])
    maj_ct <- tapply(ct_vec, cl_vec, function(v) {
      v <- v[!is.na(v) & nzchar(v)]
      if (!length(v)) NA_character_ else names(sort(table(v), decreasing = TRUE))[1]
    })

    # Backup & compute marker genes per Seurat cluster (only.pos, same thresholds)
    id_backup <- Seurat::Idents(obj)
    on.exit({try(Seurat::Idents(obj) <- id_backup, silent = TRUE)}, add = TRUE)
    Seurat::Idents(obj) <- "seurat_clusters"
    cl_markers_all <- tryCatch(
      Seurat::FindAllMarkers(
        obj, assay = assay, only.pos = TRUE,
        logfc.threshold = 0.25, min.pct = 0.2, verbose = FALSE
      ),
      error = function(e) NULL
    )

    if (is.null(cl_markers_all) || !nrow(cl_markers_all)) {
      utils::write.csv(data.frame(), cluster_go_csv, row.names = FALSE)
    } else {
      topN <- max(1L, as.integer(top_n))

      cl_top <- lapply(split(cl_markers_all, cl_markers_all$cluster), function(df) {
        lfc_col <- intersect(c("avg_log2FC","avg_logFC"), colnames(df))[1]
        if (is.na(lfc_col)) lfc_col <- setdiff(colnames(df), c("gene","cluster"))[1]
        p_col   <- if ("p_val_adj" %in% colnames(df)) "p_val_adj"
        else if ("p_val" %in% colnames(df)) "p_val" else NULL
        ord <- if (!is.null(p_col)) order(-df[[lfc_col]], df[[p_col]], na.last = TRUE)
        else                        order(-df[[lfc_col]],                na.last = TRUE)
        head(df$gene[ord], topN)
      })

      # Universe & SYMBOL->ENTREZ mapping
      all_symbols <- rownames(Seurat::GetAssayData(obj, assay = assay, slot = "data"))
      sym2ent_all <- tryCatch(
        clusterProfiler::bitr(all_symbols, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = OrgDb),
        error = function(e) NULL
      )
      universe_entrez <- if (is.null(sym2ent_all)) character(0) else unique(na.omit(sym2ent_all$ENTREZID))

      all_clusters <- as.character(sort(unique(obj$seurat_clusters)))
      iter_clusters <- unique(c(all_clusters, names(cl_top)))

      if (!exists("id_from",    inherits = FALSE)) id_from  <- "SYMBOL"
      if (!exists("strategy",   inherits = FALSE)) strategy <- "raw"
      if (!exists("key2ent",    inherits = FALSE)) key2ent  <- setNames(character(0), character(0))

      genes_union <- unique(unlist(cl_top, use.names = FALSE))
      key2ent2 <- key2ent

      # add any missing cluster genes via the same strategy-aware mapper
      missing <- setdiff(genes_union, names(key2ent2))
      if (length(missing)) {
        add <- .scq_try_map_strategy(missing, OrgDb, id_from, strategy = strategy)
        if (is.data.frame(add) && nrow(add) && "ENTREZID" %in% names(add)) {
          kcol <- intersect(c(id_from,"SYMBOL","ALIAS","ENSEMBL","ENSEMBLTRANS","ENSEMBLPROT"), colnames(add))
          kcol <- if (length(kcol)) kcol[1] else "SYMBOL"
          add_map <- setNames(add$ENTREZID, as.character(add[[kcol]]))
          add_map <- add_map[setdiff(names(add_map), names(key2ent2))]
          key2ent2 <- c(key2ent2, add_map)
        }
      }
      if (!length(key2ent2)) key2ent2 <- setNames(character(0), character(0))

      go_rows <- list()

      for (cl in iter_clusters) {
        genes <- unique(cl_top[[cl]])

        if (is.null(genes) || !length(genes)) {
          go_rows[[cl]] <- data.frame(
            cluster = cl,
            ID = NA_character_,
            Description = "(no significant BP term)",
            p.adjust = NA_real_,
            Count = 0L,
            stringsAsFactors = FALSE
          )
          if (exists(".dbg", mode="function")) .dbg(sprintf("[GO/cluster] C%s: no markers -> placeholder", cl))
          next
        }

        # map via dictionary; if empty, attempt a small robust map
        g_entrez <- unique(stats::na.omit(unname(key2ent2[as.character(genes)])))
        if (!length(g_entrez)) {
          tmp <- .scq_try_map_strategy(genes, OrgDb, id_from, strategy = strategy)
          if (is.data.frame(tmp) && nrow(tmp) && "ENTREZID" %in% names(tmp)) {
            g_entrez <- unique(stats::na.omit(tmp$ENTREZID))
          }
        }

        if (exists(".dbg", mode="function")) .dbg(sprintf("[GO/cluster] C%s: top=%d -> mapped=%d | universe=%d",
                                                          cl, length(genes), length(g_entrez), length(universe_entrez)))

        if (!length(g_entrez) || length(universe_entrez) < 100) {
          go_rows[[cl]] <- data.frame(
            cluster = cl,
            ID = NA_character_,
            Description = "(no significant BP term)",
            p.adjust = NA_real_,
            Count = 0L,
            stringsAsFactors = FALSE
          )
          if (exists(".dbg", mode="function")) .dbg(sprintf("[GO/cluster] C%s: skipped (too few IDs or small universe)", cl))
          next
        }

        ego <- tryCatch(
          suppressMessages(
            clusterProfiler::enrichGO(
              gene          = g_entrez,
              universe      = universe_entrez,
              OrgDb         = OrgDb,
              keyType       = "ENTREZID",
              ont           = "BP",
              pAdjustMethod = "BH",
              readable      = TRUE
            )
          ),
          error = function(e) { if (exists(".dbg", mode="function")) .dbg(sprintf("[GO/cluster] C%s: enrichGO ERROR: %s", cl, e$message)); NULL }
        )

        if (is.null(ego) || is.null(ego@result) || nrow(ego@result) == 0) {
          go_rows[[cl]] <- data.frame(
            cluster = cl,
            ID = NA_character_,
            Description = "(no significant BP term)",
            p.adjust = NA_real_,
            Count = 0L,
            stringsAsFactors = FALSE
          )
          if (exists(".dbg", mode="function")) .dbg(sprintf("[GO/cluster] C%s: no terms", cl))
        } else {
          df <- as.data.frame(ego)
          df <- df[order(df$p.adjust, -df$Count), , drop = FALSE]
          if (NROW(df) > 0) {
            df$cluster <- rep(cl, NROW(df))  # <-- SAFE ADD
            go_rows[[cl]] <- utils::head(df[, c("cluster","ID","Description","p.adjust","Count")], top_go_per_cluster)
          } else {
            go_rows[[cl]] <- data.frame(
              cluster = cl,
              ID = NA_character_,
              Description = "(no significant BP term)",
              p.adjust = NA_real_,
              Count = 0L,
              stringsAsFactors = FALSE
            )
          }
        }
      }

      go_cl <- do.call(rbind, go_rows)
      if (!is.null(go_cl) && nrow(go_cl)) {
        go_cl <- as.data.frame(go_cl, stringsAsFactors = FALSE)
        names(go_cl)[names(go_cl) == "Description"] <- "term"
        go_cl$p_neglog10 <- -log10(go_cl$p.adjust)
        cl_order <- sort(unique(as.integer(as.character(go_cl$cluster))))
        go_cl$cluster <- factor(as.character(go_cl$cluster), levels = as.character(cl_order))
        lab_celltype <- maj_ct[as.character(go_cl$cluster)]
        lab_celltype[is.na(lab_celltype) | !nzchar(lab_celltype)] <- "Unknown"
        go_cl$term_label <- paste0("C", as.character(go_cl$cluster), ": ", lab_celltype, " _ ", go_cl$term)
        ord <- order(as.integer(as.character(go_cl$cluster)), -go_cl$p_neglog10, -go_cl$Count)
        go_cl <- go_cl[ord, , drop = FALSE]
        go_cl$term_label <- factor(go_cl$term_label, levels = rev(unique(go_cl$term_label)))

        utils::write.csv(go_cl, cluster_go_csv, row.names = FALSE)

        if (requireNamespace("ggplot2", quietly = TRUE)) {
          gg <- ggplot2::ggplot(
            go_cl,
            ggplot2::aes(x = p_neglog10, y = term_label, fill = cluster)
          ) +
            ggplot2::geom_col() +
            ggplot2::labs(
              x = expression(-log[10]("adj. P")),
              y = NULL, fill = "Cluster",
              title = "Top GO BP term per Seurat cluster"
            ) +
            ggplot2::theme_minimal(base_size = 10) +
            ggplot2::theme(
              legend.position = "none",
              panel.grid = ggplot2::element_blank(),
              axis.text.y = ggplot2::element_text(size = 9)
            )
          ggplot2::ggsave(
            cluster_go_png, gg,
            width = 8.5, height = nrow(go_cl) * 0.15 + 1, dpi = 150, bg = "transparent"
          )
        }
      } else {
        utils::write.csv(data.frame(), cluster_go_csv, row.names = FALSE)
      }
    }
  } else {
    utils::write.csv(data.frame(), cluster_go_csv, row.names = FALSE)
  }

  invisible(list(heatmap = heat_png, enrich = go_png))
}







validation_plots2_post_annotation <- function(
    obj,
    save_dir       = "qc_outputs",
    assay          = Seurat::DefaultAssay(obj),
    species        = c("mouse","human"),
    panglao_file   = system.file("extdata","PanglaoDB_markers_27_Mar_2020.tsv",
                                 package = utils::packageName()),
    panglao_ui_max = 0.20,
    fm_top_n       = 50,      # DE genes considered per cluster *before* filter
    fm_padj_max    = 0.05,    # padj cutoff for DE rows
    per_cluster_keep = 10,    # genes kept *after* Panglao filter for the plot
    overlap_mode   = c("count","jaccard"),
    width_in       = 9, height_in = 7, dpi = 150
){
  stopifnot(inherits(obj,"Seurat"))
  species <- match.arg(species); overlap_mode <- match.arg(overlap_mode)

  out_png <- file.path(save_dir, "celltype_gene_heatmap2.png")
  dbgdir  <- file.path(save_dir, "annot_validation")
  dir.create(dbgdir, recursive = TRUE, showWarnings = FALSE)

  .log <- function(fmt, ...) {
    ln <- if (length(list(...))) do.call(sprintf, c(list(fmt), list(...))) else fmt
    cat(sprintf("[%s] %s\n", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), ln),
        file = file.path(dbgdir, "log.txt"), append = TRUE)
    if (isTRUE(getOption("scQCenrich.debug", FALSE))) message(ln)
  }
  .csv <- function(name, x, rn = TRUE) { p <- file.path(dbgdir, name); try(utils::write.csv(x, p, row.names = rn), silent=TRUE); p }
  .snippet <- function(lines) writeLines(lines, file.path(dbgdir, "RESULT_SNIPPET.txt"))

  message("[start] validation_plots2_post_annotation | assay=%s species=%s", assay, species)
  if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)

  # --- 0) Ensure normalized layer exists (no mutation beyond 'data' if missing)
  Seurat::DefaultAssay(obj) <- assay
  lyr <- try(Seurat::GetAssayData(obj, assay = assay, layer = "data"), silent = TRUE)
  if (inherits(lyr, "try-error") || !nrow(lyr) || !ncol(lyr)) {
    message("[prep] NormalizeData(): creating 'data' layer")
    obj <- Seurat::NormalizeData(obj, verbose = FALSE)
  }

  # --- 1) Load Panglao *canonical* signatures (UI filter)
  if (!nzchar(panglao_file) || !file.exists(panglao_file)) {
    message("[abort] Panglao file not found: %s", panglao_file)
  }
  sig0 <- try(panglao_signatures(tsv = panglao_file, species = species,
                                 canonical_only = TRUE, ui_max = panglao_ui_max,
                                 min_genes = 2), silent = TRUE)
  if (inherits(sig0, "try-error") || !length(sig0)) {
    message("[warn] 0 Panglao signatures available after filters; will still build a DE-only heatmap.")
    sig0 <- list()
  } else {
    message("[panglao] signatures kept: %d (canonical=TRUE ui_max=%.2f)", length(sig0), panglao_ui_max)
  }

  # Map markers to object feature case
  genes_obj <- rownames(Seurat::GetAssayData(obj, assay = assay, layer = "data"))
  if (!length(genes_obj)) genes_obj <- rownames(Seurat::GetAssayData(obj, assay = assay, layer = "counts"))
  up_obj <- toupper(genes_obj)
  map_to_obj <- function(v){
    idx <- match(toupper(v), up_obj); idx <- idx[!is.na(idx)]
    if (!length(idx)) character(0) else genes_obj[idx]
  }
  sig_obj <- lapply(sig0, map_to_obj)
  cov_df <- data.frame(
    celltype = names(sig0),
    n_markers = vapply(sig0, length, 1L),
    overlap   = vapply(sig_obj, length, 1L),
    frac = vapply(sig_obj, function(x) ifelse(length(x), length(x)/max(1, length(x)), 0), 1),
    check.names = FALSE
  )
  .csv("01_panglao_coverage.csv", cov_df)

  # --- 2) Ensure clusters, then DE per cluster (pos only)
  if (!("seurat_clusters" %in% colnames(obj@meta.data))) {
    message("[prep] computing seurat_clusters for DE")
    obj <- Seurat::FindVariableFeatures(obj, nfeatures = 2000, verbose = FALSE)
    obj <- Seurat::ScaleData(obj, verbose = FALSE)
    obj <- Seurat::RunPCA(obj, npcs = 30, verbose = FALSE)
    obj <- Seurat::FindNeighbors(obj, dims = 1:20, verbose = FALSE)
    obj <- Seurat::FindClusters(obj, resolution = 0.4, verbose = FALSE)
  } else {
    message("[prep] using existing seurat_clusters (K=%d)", length(unique(obj$seurat_clusters)))
  }
  Seurat::Idents(obj) <- "seurat_clusters"

  de <- try(Seurat::FindAllMarkers(obj, assay = assay, only.pos = TRUE,
                                   logfc.threshold = 0.25, min.pct = 0.10, verbose = FALSE),
            silent = TRUE)
  if (inherits(de, "try-error") || !is.data.frame(de) || !nrow(de)) {
    message("[abort] FindAllMarkers yielded 0 rows; writing placeholder image.")
    .placeholder_png(out_png, "No DE markers found.\nCannot build validation heatmap.")
    return(invisible(out_png))
  }
  fc_col   <- intersect(c("avg_log2FC","avg_logFC"), colnames(de))[1]
  padj_col <- intersect(c("p_val_adj","p_val"), colnames(de))[1]
  de <- de[order(de$cluster, -abs(de[[fc_col]]), de[[padj_col]]), , drop = FALSE]
  if (!is.na(padj_col)) de <- de[is.na(de[[padj_col]]) | de[[padj_col]] <= fm_padj_max, , drop = FALSE]
  .csv("02_findallmarkers_head.csv", utils::head(de, 50))
  message("[DE] rows after filters: %d (topN per cluster=%d, padj<=%.3g)", nrow(de), fm_top_n, fm_padj_max)

  # --- 3) Decide "respective cell type" per cluster
  cl_vec <- as.character(obj$seurat_clusters); names(cl_vec) <- colnames(obj)
  clusters <- sort(unique(cl_vec))

  # First choice: majority of per-cell labels if present
  ctype_by_cluster <- setNames(rep(NA_character_, length(clusters)), clusters)
  if ("auto_celltype" %in% colnames(obj@meta.data)) {
    tab <- tapply(obj$auto_celltype, cl_vec, function(v) {
      v <- as.character(v); v <- v[nzchar(v) & !is.na(v)]
      if (!length(v)) NA_character_ else names(sort(table(v), decreasing = TRUE))[1]
    })
    ctype_by_cluster[names(tab)] <- as.character(tab)
  }

  # Fallback: best overlap between top DE genes and Panglao signatures
  U <- function(x) unique(toupper(x))
  de_by_cl <- split(as.character(de$gene), de$cluster)
  de_by_cl <- lapply(de_by_cl, function(v) unique(utils::head(v, fm_top_n)))
  if (length(sig_obj)) {
    for (k in names(ctype_by_cluster)) {
      if (!is.na(ctype_by_cluster[[k]])) next
      topU <- U(de_by_cl[[k]])
      sc <- sapply(sig_obj, function(ref) {
        inter <- length(intersect(topU, U(ref)))
        if (identical(overlap_mode, "jaccard"))
          inter/length(union(topU, U(ref))) else inter
      })
      best <- names(sc)[which.max(sc)]
      if (length(sc) && is.finite(max(sc)) && max(sc) > 0) ctype_by_cluster[[k]] <- best
    }
  }
  map_df <- data.frame(cluster = names(ctype_by_cluster), celltype = ctype_by_cluster, check.names = FALSE)
  .csv("03_cluster_celltype_map.csv", map_df)

  # --- 4) For each cluster, keep only DE genes that are canonical markers of its celltype
  keep_by_cl <- list()
  kept_counts <- integer(length(clusters)); names(kept_counts) <- clusters
  for (k in clusters) {
    top_genes <- de_by_cl[[k]] %||% character(0)
    ct <- ctype_by_cluster[[as.character(k)]]
    if (!is.na(ct) && ct %in% names(sig_obj)) {
      g <- intersect(map_to_obj(top_genes), sig_obj[[ct]])
      g <- unique(utils::head(g, per_cluster_keep))
    } else {
      g <- character(0)
    }
    # graceful fallback: if nothing left, keep top 3 DE to avoid blank figure
    if (!length(g)) g <- unique(utils::head(map_to_obj(top_genes), min(3, per_cluster_keep)))
    keep_by_cl[[as.character(k)]] <- g
    kept_counts[[as.character(k)]] <- length(g)
  }
  kept_tbl <- data.frame(cluster = names(keep_by_cl),
                         n_genes = vapply(keep_by_cl, length, 1L),
                         celltype = ctype_by_cluster[names(keep_by_cl)],
                         check.names = FALSE)
  .csv("04_genes_kept_per_cluster.csv", kept_tbl)

  # Union gene set & build cluster-averaged expression matrix
  genes_keep <- unique(unlist(keep_by_cl, use.names = FALSE))
  if (!length(genes_keep)) {
    message("[warn] After filtering, 0 genes remained; writing placeholder PNG.")
    .placeholder_png(out_png, "No overlap with Panglao canonical markers.\n(See annot_validation/log.txt)")
    return(invisible(out_png))
  }
  M <- Seurat::GetAssayData(obj, assay = assay, layer = "data")
  M <- M[intersect(genes_keep, rownames(M)), , drop = FALSE]
  if (!nrow(M)) {
    .placeholder_png(out_png, "Selected genes not present in object.")
    return(invisible(out_png))
  }
  avg_by_cl <- sapply(clusters, function(k) Matrix::rowMeans(M[, cl_vec == k, drop = FALSE]))
  colnames(avg_by_cl) <- paste0("C", clusters, ifelse(is.na(ctype_by_cluster), "", paste0(":", ctype_by_cluster)))
  # z-score per gene for prettier heatmap
  z <- t(scale(t(as.matrix(avg_by_cl))))
  z[!is.finite(z)] <- 0
  .csv("05_heatmap_matrix_head.csv", utils::head(as.data.frame(z), 20), rn = TRUE)

  # --- 5) Plot PNG (ComplexHeatmap if available; else ggplot)
  if (requireNamespace("ComplexHeatmap", quietly = TRUE)) {
    hm <- ComplexHeatmap::Heatmap(z,
                                  name = "z",
                                  column_title = "Top DE markers filtered by Panglao canonical markers",
                                  row_names_gp = grid::gpar(cex = 0.6),
                                  column_names_gp = grid::gpar(cex = 0.7))
    grDevices::png(out_png, width = width_in, height = height_in, units = "in", res = dpi, bg = getOption("scQCenrich.plot_bg","white"))
    print(hm); grDevices::dev.off()
  } else {
    df <- as.data.frame(z)
    df$gene <- rownames(df)
    long <- reshape2::melt(df, id.vars = "gene", variable.name = "cluster", value.name = "z")
    p <- ggplot2::ggplot(long, ggplot2::aes(x = cluster, y = gene, fill = z)) +
      ggplot2::geom_tile() +
      ggplot2::scale_fill_viridis_c() +
      ggplot2::labs(title = "Top DE markers filtered by Panglao canonical markers", x = NULL, y = NULL) +
      ggplot2::theme_minimal(base_size = 10) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
    ggplot2::ggsave(out_png, p, width = width_in, height = height_in, dpi = dpi,
                    bg = getOption("scQCenrich.plot_bg","white"))
  }
  message("[done] wrote %s (genes=%d, clusters=%d)", out_png, nrow(z), ncol(z))

  # --- 6) Result snippet (compact, easy to read in CI logs)
  preview <- paste(utils::head(rownames(z), 8), collapse = ", ")
  labline <- paste(utils::head(colnames(z), 6), collapse = " | ")
  unk_rate <- mean(is.na(ctype_by_cluster)) * 100
  .snippet(c(
    "==== Annotation validation (markers) ====",
    sprintf("clusters: %d | genes plotted: %d | unk celltype rate: %.1f%%", length(clusters), nrow(z), unk_rate),
    sprintf("padj<=%.3g, topN=%d, per_cluster_keep=%d, mode=%s", fm_padj_max, fm_top_n, per_cluster_keep, overlap_mode),
    sprintf("clusters(head): %s", labline),
    sprintf("genes(head): %s", preview),
    "========================================"
  ))

  invisible(out_png)
}

.placeholder_png <- function(path, msg) {
  grDevices::png(path, width = 1000, height = 700, res = 150, bg = getOption("scQCenrich.plot_bg","white"))
  op <- par(mar = c(0,0,0,0)); on.exit(par(op), add = TRUE)
  plot.new(); text(0.5, 0.5, msg, cex = 1.1); grDevices::dev.off()
}




#' FeaturePlot-like map of auto-annotation confidence
#' @param obj Seurat object
#' @param imgDir directory to write the PNG (created if missing)
#' @param reduction which embedding to use; default: first of umap/tsne/pca found
#' @param col meta column name for confidence; default tries common names
#' @param filename output filename
#' @return full path to the PNG (character) or NULL if skipped
#' @export
scq_plot_auto_confidence <- function(
    obj,
    imgDir   = "qc_outputs",
    reduction = NULL,
    col       = NULL,
    filename  = "auto_confidence.png"
){
  stopifnot(inherits(obj, "Seurat"))
  if (!dir.exists(imgDir)) dir.create(imgDir, recursive = TRUE, showWarnings = FALSE)

  # choose reduction
  if (is.null(reduction)) {
    reduction <- intersect(c("umap","tsne","pca"), names(obj@reductions))
    reduction <- if (length(reduction)) reduction[[1]] else NULL
  }
  if (is.null(reduction)) {
    message("[scQCenrich] scq_plot_auto_confidence(): no embedding found; skipping.")
    return(NULL)
  }

  # choose confidence column
  if (is.null(col)) {
    candidates <- c("auto_confidence","auto_score","sctype_confidence","sctype_score",
                    "singleR_score","label_transfer_score","auto_celltype_conf")
    col <- candidates[candidates %in% colnames(obj@meta.data)]
    col <- if (length(col)) col[[1]] else NULL
  }
  if (is.null(col)) {
    message("[scQCenrich] scq_plot_auto_confidence(): no confidence column found; skipping.")
    return(NULL)
  }

  vals <- obj@meta.data[[col]]
  if (!is.numeric(vals)) {
    message("[scQCenrich] scq_plot_auto_confidence(): confidence column is not numeric; skipping.")
    return(NULL)
  }

  emb <- Seurat::Embeddings(obj[[reduction]])
  df  <- data.frame(x = emb[,1], y = emb[,2], conf = vals, stringsAsFactors = FALSE)

  p <- ggplot2::ggplot(df, ggplot2::aes(x = x, y = y, color = conf)) +
    ggplot2::geom_point(size = 0.6, alpha = 0.8) +
    ggplot2::labs(title = sprintf("Auto-annotation confidence (%s)", reduction),
                  color = "confidence", x = reduction %||% "dim1", y = reduction %||% "dim2") +
    ggplot2::theme_minimal()

  outfile <- file.path(imgDir, filename)

  # robust PNG writer (uses scq_png if present, else base png)
  if (exists("scq_png", mode = "function")) {
    scq_png(outfile, width = 1200, height = 900, res = 120, expr = {
      print(p)
    })
  } else {
    grDevices::png(outfile, width = 1200, height = 900, res = 120)
    print(p)
    grDevices::dev.off()
  }

  return(outfile)
}

`%||%` <- function(a, b) if (!is.null(a)) a else b
