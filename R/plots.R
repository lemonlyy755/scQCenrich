# ---- plots.R ----
qc_diagnostic_panel <- function(
    obj, metrics, status_df,
    assay = DefaultAssay(obj),
    save_dir = "qc_outputs",
    prefix = "qc",
    sample_col = NULL,
    group_col = "auto_celltype"  # NEW: functional grouping for bars
){
  if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)
  obj <- .ensure_pca_umap(obj, assay = assay, npcs = 20)

  # --- 1) UMAP by QC status ---
  umap_df <- .make_umap_df(obj, color_col = "qc_status")
  p_qc <- ggplot2::ggplot(umap_df, ggplot2::aes(UMAP_1, UMAP_2, color = qc_status)) +
    ggplot2::geom_point(size = 0.3, alpha = 0.6) + ggplot2::theme_classic() +
    ggplot2::labs(title = "UMAP (colored by QC status)") +
    ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = 3, alpha = 1)))
  f_qc <- file.path(save_dir, sprintf("%s_umap_qc.png", prefix)); .gg_save_safe(p_qc, f_qc, 6, 5, 200)

  # --- 1b) UMAP by annotation (if present) ---
  annot_col <- if (!is.null(group_col) && group_col %in% colnames(obj@meta.data)) group_col else NULL
  f_annot <- NA_character_
  if (!is.null(annot_col)) {
    ua <- .make_umap_df(obj, color_col = annot_col)
    p_annot <- ggplot2::ggplot(ua, ggplot2::aes(UMAP_1, UMAP_2, color = .data[[annot_col]])) +
      ggplot2::geom_point(size = 0.3, alpha = 0.6) + ggplot2::theme_classic() +
      ggplot2::labs(title = paste0("UMAP (", annot_col, ")"), color = annot_col) +
      ggplot2::guides(color = ggplot2::guide_legend(override.aes = list(size = 3, alpha = 1)))
    f_annot <- file.path(save_dir, sprintf("%s_umap_annotation.png", prefix)); .gg_save_safe(p_annot, f_annot, 7, 5, 200)
  }

  # --- 2) UMAP by doublet ---
  md <- obj@meta.data; md$cell <- rownames(md)
  dd <- .make_umap_df(obj); dd$is_doublet <- md[rownames(dd), "is_doublet", drop = TRUE]
  p_dbl <- ggplot2::ggplot(dd, ggplot2::aes(UMAP_1, UMAP_2, color = is_doublet)) +
    ggplot2::geom_point(size = 0.3, alpha = 0.6) + ggplot2::theme_classic() +
    ggplot2::scale_color_manual(values = c("FALSE"="#A0A0A0","TRUE"="#E74C3C"), na.translate = TRUE) +
    ggplot2::labs(title = "UMAP (predicted doublets)")
  f_dbl <- file.path(save_dir, sprintf("%s_umap_doublet.png", prefix)); .gg_save_safe(p_dbl, f_dbl, 6, 5, 200)

  # --- 3) Doublet rate bar chart by functional group (robust) ---
  md <- obj@meta.data

  # Coerce doublet flag robustly (handles TRUE/FALSE, 0/1, "true"/"false", etc.)
  is_dbl <- md[["is_doublet"]]
  if (is.null(is_dbl)) is_dbl <- NA
  if (!is.logical(is_dbl)) {
    if (is.numeric(is_dbl)) {
      is_dbl <- is_dbl > 0
    } else {
      is_dbl <- tolower(as.character(is_dbl)) %in% c("true","t","1","yes")
    }
  }
  md$.__is_doublet__ <- is_dbl

  gcol <- if (!is.null(annot_col)) {
    annot_col
  } else if (!is.null(sample_col) && sample_col %in% colnames(md)) {
    sample_col
  } else {
    NA_character_
  }

  if (is.na(gcol)) {
    rate <- mean(md$.__is_doublet__, na.rm = TRUE)
    tab  <- data.frame(group = "All cells",
                       rate  = ifelse(is.finite(rate), rate, NA_real_))
  } else {
    md$.__group__ <- md[[gcol]]
    md2 <- md[!is.na(md$.__group__), , drop = FALSE]

    if (nrow(md2)) {
      # use 'by=' form to avoid formula NSE headaches
      tab <- aggregate(x  = list(rate = as.numeric(md2$.__is_doublet__)),
                       by = list(group = as.character(md2$.__group__)),
                       FUN = function(x) mean(x == 1, na.rm = TRUE))
    } else {
      tab <- data.frame(group = character(0), rate = numeric(0))
    }

    # Fallback if everything was NA
    if (!nrow(tab)) {
      rate <- mean(md$.__is_doublet__, na.rm = TRUE)
      tab  <- data.frame(group = "All cells",
                         rate  = ifelse(is.finite(rate), rate, NA_real_))
    }
  }

  tab <- tab[order(tab$rate, decreasing = TRUE), , drop = FALSE]

  p_bar <- ggplot2::ggplot(tab, ggplot2::aes(x = group, y = 100*rate)) +
    ggplot2::geom_col() + ggplot2::theme_classic() +
    ggplot2::ylab("% predicted doublets") +
    ggplot2::xlab(if (is.na(gcol)) NULL else gcol) +
    ggplot2::coord_cartesian(ylim = c(0, max(5, 100*max(tab$rate, na.rm = TRUE)))) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1))
  f_bar <- file.path(save_dir, sprintf("%s_doublet_rate_by_group.png", prefix))
  .gg_save_safe(p_bar, f_bar, 8.5, 4.5, 200)

  invisible(list(umap_qc = f_qc, umap_annot = f_annot, umap_doublet = f_dbl, bar_doublet = f_bar))
}


.ensure_pca_umap <- function(obj, assay = DefaultAssay(obj), npcs = 20){
  odef <- DefaultAssay(obj); on.exit(DefaultAssay(obj) <- odef, add = TRUE)
  DefaultAssay(obj) <- assay

  # If 'data' layer is empty, create it for diagnostics
  lyr_data <- try(Seurat::GetAssayData(obj, assay = assay, layer = "data"), silent = TRUE)
  if (inherits(lyr_data, "try-error") || nrow(lyr_data) == 0 || ncol(lyr_data) == 0) {
    if (isTRUE(getOption("scQCenrich.debug"))) message("[scQCenrich-debug] NormalizeData(): creating non-empty 'data' layer for diagnostics. ")
    obj <- Seurat::NormalizeData(obj, verbose = FALSE)
  }

  if (!("pca" %in% names(obj@reductions))) {
    if (isTRUE(getOption("scQCenrich.debug"))) message("[scQCenrich-debug] FindVariableFeatures()/ScaleData()/RunPCA(): creating PCA. ")
    obj <- Seurat::FindVariableFeatures(obj, nfeatures = 2000, verbose = FALSE)
    obj <- Seurat::ScaleData(obj, verbose = FALSE)
    obj <- Seurat::RunPCA(obj, npcs = npcs, verbose = FALSE)
  }
  if (!("umap" %in% names(obj@reductions)) || ncol(Seurat::Embeddings(obj, "umap")) < 2) {
    if (isTRUE(getOption("scQCenrich.debug"))) message("[scQCenrich-debug] RunUMAP(): creating UMAP with 20 PCs. ")
    nd <- min(npcs, ncol(Seurat::Embeddings(obj, "pca")))
    obj <- Seurat::RunUMAP(
      obj, reduction = "pca", dims = 1:nd, verbose = FALSE,
      seed.use = .umap_seed()
    )

  }
  obj
}

.make_umap_df <- function(obj, color_col = NULL){
  emb <- Seurat::Embeddings(obj, "umap")
  stopifnot(ncol(emb) >= 2)
  df <- data.frame(
    cell    = rownames(emb),
    UMAP_1  = emb[, 1],
    UMAP_2  = emb[, 2],
    row.names = rownames(emb),
    check.names = FALSE
  )
  if (!is.null(color_col) && color_col %in% colnames(obj@meta.data)) {
    df[[color_col]] <- obj@meta.data[rownames(df), color_col, drop = TRUE]
  }
  df
}

.gg_save_safe <- function(plot_or_path, filename, width = 7, height = 5, dpi = 150,
                          bg = getOption("scQCenrich.plot_bg", "white")) {
  if (inherits(plot_or_path, "gg")) {
    ggplot2::ggsave(filename, plot = plot_or_path, width = width, height = height, dpi = dpi, bg = bg)
  } else if (is.character(plot_or_path) && file.exists(plot_or_path)) {
    return(invisible(filename))
  }
  invisible(filename)
}


# ---- plots.R (append) ----

.winsorize <- function(x, q = c(0.02, 0.98)) {
  lo <- stats::quantile(x, q[1], na.rm = TRUE)
  hi <- stats::quantile(x, q[2], na.rm = TRUE)
  x[x < lo] <- lo; x[x > hi] <- hi
  x
}

.pretty_feature_plot <- function(obj, value, title, point_size = 0.25, alpha = 0.7) {
  stopifnot(length(value) == ncol(obj))
  df <- .make_umap_df(obj)
  df$val <- value[rownames(df)]
  df <- df[is.finite(df$val), , drop = FALSE]

  BG <- getOption("scQCenrich.plot_bg", "white")

  ggplot2::ggplot(df, ggplot2::aes(UMAP_1, UMAP_2, color = val)) +
    ggplot2::geom_point(size = point_size, alpha = alpha) +
    ggplot2::theme_void(base_size = 10) +
    ggplot2::theme(
      plot.background   = ggplot2::element_rect(fill = BG, colour = NA),
      panel.background  = ggplot2::element_rect(fill = BG, colour = NA),
      legend.background = ggplot2::element_rect(fill = BG, colour = NA),
      legend.key        = ggplot2::element_rect(fill = BG, colour = NA),
      legend.position   = "right",              # <-- per-plot legend on the right
      legend.margin     = ggplot2::margin(2,2,2,2,"mm"),
      legend.box.margin = ggplot2::margin(2,4,2,2,"mm"),
      plot.title        = ggplot2::element_text(hjust = 0.5, face = "bold"),
      plot.margin       = ggplot2::margin(2,6,2,2,"mm")
    ) +
    ggplot2::labs(title = title, color = NULL) +
    ggplot2::scale_color_viridis_c() +
    ggplot2::guides(color = ggplot2::guide_colorbar(
      title.position = "top",
      barheight = grid::unit(16, "mm"),
      barwidth  = grid::unit(3,  "mm")
    ))
}



#' Make & save UMAP featuremaps for QC metrics
#' @param obj Seurat object (must have UMAP). Use .ensure_pca_umap() beforehand.
#' @param metrics Vector of meta column names to plot (only those present are used)
#' @param save_dir Where to write PNGs
#' @param prefix File prefix for images
#' @param ncol Number of columns in the grid image
#' @param logify Named logical vector: which metrics to log1p/ log10; defaults below
#' @export
qc_featuremaps <- function(
    obj,
    metrics = c("nCount_QC","nFeature_QC","pctMT_QC","MALAT1_frac_QC",
                "stress_score_QC","intronic_frac_QC",
                "spliced_frac_QC","unspliced_frac_QC","u2s_ratio_QC","dbl_score"),
    save_dir = "qc_outputs",
    prefix   = "qc",
    ncol     = 4,
    logify   = c(nCount_QC=TRUE, nFeature_QC=TRUE, u2s_ratio_QC=TRUE, dbl_score=FALSE)
){
  if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)
  obj <- .ensure_pca_umap(obj)
  BG <- getOption("scQCenrich.plot_bg", "white")

  present <- intersect(metrics, colnames(obj@meta.data))
  if (!length(present)) return(invisible(NULL))

  # nice display labels
  labmap <- c(
    nCount_QC        = "Total counts (log)",
    nFeature_QC      = "Detected genes (log)",
    pctMT_QC         = "% mito",
    MALAT1_frac_QC   = "MALAT1 fraction",
    stress_score_QC  = "Stress score",
    intronic_frac_QC = "Intronic fraction",
    spliced_frac_QC  = "Spliced fraction",
    unspliced_frac_QC= "Unspliced fraction",
    u2s_ratio_QC     = "Unspliced/Spliced (log)",
    dbl_score        = "Doublet score"
  )

  plots <- list()
  for (m in present) {
    v <- obj@meta.data[, m, drop = TRUE]
    # transform & clip
    if (isTRUE(logify[m])) {   # <- safe: missing names become NA, isTRUE(NA) == FALSE
      if (all(v >= 0, na.rm = TRUE)) v <- log10(v + 1) else v <- log1p(v)
    }
    v <- .winsorize(v, c(0.02, 0.98))
    p <- .pretty_feature_plot(obj, setNames(v, rownames(obj@meta.data)), labmap[[m]] %||% m)
    f <- file.path(save_dir, sprintf("%s_featuremap_%s.png", prefix, m))
    .gg_save_safe(p, f, width = 6, height = 5, dpi = 220, bg = BG)
    plots[[m]] <- p
  }

  # grid image (patchwork)
  if (length(plots)) {
    if (!requireNamespace("patchwork", quietly = TRUE)) {
      warning("Install 'patchwork' to save a featuremap grid; individual PNGs were saved.")
    } else {
      nrow <- ceiling(length(plots) / ncol)

      # IMPORTANT: keep per-plot guides (no collection)
      grid <- patchwork::wrap_plots(plots, ncol = ncol, byrow = TRUE, guides = "keep") +
        patchwork::plot_annotation(
          theme = ggplot2::theme(
            plot.background  = ggplot2::element_rect(fill = BG, colour = NA),
            panel.background = ggplot2::element_rect(fill = BG, colour = NA)
          )
        )

      fgrid <- file.path(save_dir, sprintf("%s_featuremaps_grid.png", prefix))
      .gg_save_safe(grid, fgrid, width = 7.2 * ncol, height = 5 * nrow, dpi = 200, bg = BG)
    }
  }

  invisible(invisible(NULL))
}
