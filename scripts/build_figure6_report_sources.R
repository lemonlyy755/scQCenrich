#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(ggplot2)
  library(dplyr)
  library(clusterProfiler)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(patchwork)
})

options(stringsAsFactors = FALSE, digits = 17, scipen = 999)

repo_dir <- "C:/RItems/Package/scQCenrich"
out_dir <- "C:/Users/ken/Documents/Package paper/report/source"
verify_dir <- file.path(out_dir, "verification")
component_dir <- file.path(out_dir, "Figure6_report_components")
work_dir <- file.path(out_dir, "Figure6_regeneration_work")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(verify_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(component_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(work_dir, recursive = TRUE, showWarnings = FALSE)
setwd(work_dir)

if (requireNamespace("pkgload", quietly = TRUE)) {
  pkgload::load_all(repo_dir, quiet = TRUE)
} else {
  library(scQCenrich)
}

write_csv_base <- function(x, path) {
  utils::write.csv(x, path, row.names = FALSE, fileEncoding = "UTF-8")
  path
}

parse_ratio_decimal <- function(x) {
  vapply(strsplit(as.character(x), "/", fixed = TRUE), function(z) {
    if (length(z) != 2L) {
      return(NA_real_)
    }
    as.numeric(z[[1]]) / as.numeric(z[[2]])
  }, numeric(1))
}

winsorize <- function(x, q = c(0.02, 0.98)) {
  lo <- stats::quantile(x, q[1], na.rm = TRUE)
  hi <- stats::quantile(x, q[2], na.rm = TRUE)
  x[x < lo] <- lo
  x[x > hi] <- hi
  x
}

get_umap <- function(obj) {
  emb <- Seurat::Embeddings(obj, "umap")[, 1:2, drop = FALSE]
  colnames(emb) <- c("UMAP_1", "UMAP_2")
  data.frame(
    cell = rownames(emb),
    UMAP_1 = emb[, 1],
    UMAP_2 = emb[, 2],
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
}

feature_plot_from_source <- function(df, title) {
  ggplot(df, aes(UMAP_1, UMAP_2, color = plot_value)) +
    geom_point(size = 0.24, alpha = 0.7) +
    coord_equal() +
    theme_void(base_size = 10) +
    labs(title = title, color = NULL) +
    scale_color_viridis_c() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "right",
      plot.margin = margin(2, 6, 2, 2, "mm")
    )
}

make_heatmap_source <- function(obj, assay = "RNA", group.by = "auto_celltype",
                                slot.use = "data") {
  Seurat::Idents(obj) <- obj@meta.data[[group.by]]
  markers_all <- Seurat::FindAllMarkers(
    obj,
    assay = assay,
    slot = slot.use,
    only.pos = TRUE,
    min.pct = 0.1,
    logfc.threshold = 0.25,
    verbose = FALSE
  )
  if (!nrow(markers_all)) {
    stop("No cell-type markers found for Figure 6C/D.", call. = FALSE)
  }

  lfcol <- if ("avg_log2FC" %in% colnames(markers_all)) {
    "avg_log2FC"
  } else if ("avg_logFC" %in% colnames(markers_all)) {
    "avg_logFC"
  } else {
    "p_val_adj"
  }

  top_per_cluster10 <- dplyr::group_by(markers_all, cluster) |>
    dplyr::slice_max(order_by = .data[[lfcol]], n = 10, with_ties = FALSE) |>
    dplyr::ungroup()

  genes_use <- unique(top_per_cluster10$gene)
  expr_avg <- try(
    Seurat::AverageExpression(
      obj,
      assays = assay,
      slot = slot.use,
      features = genes_use,
      return.seurat = FALSE
    ),
    silent = TRUE
  )
  if (!inherits(expr_avg, "try-error") && is.list(expr_avg) && length(expr_avg) >= 1) {
    mat <- expr_avg[[assay]]
  } else {
    agg <- Seurat::AggregateExpression(
      obj,
      assays = assay,
      slot = slot.use,
      features = genes_use
    )
    mat <- agg[[assay]]
  }

  mat <- as.matrix(mat[intersect(rownames(mat), genes_use), , drop = FALSE])
  if (!is.null(rownames(mat))) {
    rownames(mat) <- make.unique(as.character(rownames(mat)))
  }

  norm_label <- function(x) gsub("-", "_", x, fixed = TRUE)
  colnames(mat) <- norm_label(colnames(mat))
  tpc_cluster <- norm_label(top_per_cluster10$cluster)

  md <- obj@meta.data
  all_groups <- colnames(mat)
  cell_counts_tab <- table(factor(norm_label(md[[group.by]]), levels = all_groups))
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
      ordered_genes <- c(ordered_genes, head(rownames(mat)[ord], 5))
    }
  }
  ordered_genes <- unique(ordered_genes)
  if (!length(ordered_genes)) {
    ord_global <- order(rowMeans(abs(mat), na.rm = TRUE), decreasing = TRUE)
    ordered_genes <- head(rownames(mat)[ord_global], max(1, min(30, nrow(mat))))
  }
  mat <- mat[ordered_genes, , drop = FALSE]

  src <- as.data.frame(as.table(mat), stringsAsFactors = FALSE)
  names(src) <- c("gene", "auto_celltype", "average_expression")
  src$gene_order <- match(src$gene, rownames(mat))
  src$celltype_order <- match(src$auto_celltype, colnames(mat))
  src <- src[order(src$gene_order, src$celltype_order), , drop = FALSE]

  list(
    markers_all = markers_all,
    lfcol = lfcol,
    heatmap_matrix = mat,
    heatmap_source = src,
    celltype_order = colnames(mat)
  )
}

go_sources <- function(obj, heat, assay = "RNA", group.by = "auto_celltype",
                       top_n = 100, min_genes_for_enrich = 10,
                       top_go_per_cluster = 2) {
  scq_try_map_strategy <- get(".scq_try_map_strategy", envir = asNamespace("scQCenrich"))
  markers_all <- heat$markers_all
  lfcol <- heat$lfcol

  top_per_cluster <- dplyr::group_by(markers_all, cluster) |>
    dplyr::slice_max(order_by = .data[[lfcol]], n = top_n, with_ties = FALSE) |>
    dplyr::ungroup()

  feats <- rownames(Seurat::GetAssayData(obj, assay = assay, layer = "counts"))
  kt_all <- try(AnnotationDbi::keytypes(org.Hs.eg.db::org.Hs.eg.db), silent = TRUE)
  if (inherits(kt_all, "try-error")) {
    kt_all <- character(0)
  }
  kt_cand <- intersect(
    c("ENTREZID", "SYMBOL", "ALIAS", "ENSEMBL", "ENSEMBLTRANS",
      "ENSEMBLPROT", "REFSEQ", "UNIPROT", "GENENAME"),
    kt_all
  )
  if (!length(kt_cand)) {
    kt_cand <- "SYMBOL"
  }

  probe_rows <- list()
  for (kt in kt_cand) {
    for (st in c("raw", "title", "upper", "auto")) {
      res <- try(scq_try_map_strategy(feats, org.Hs.eg.db::org.Hs.eg.db, kt, strategy = st),
        silent = TRUE
      )
      nm <- if (inherits(res, "try-error") || !nrow(res)) {
        0L
      } else {
        length(unique(stats::na.omit(res$ENTREZID)))
      }
      probe_rows[[paste(kt, st, sep = ":")]] <- data.frame(
        KEYTYPE = kt,
        STRATEGY = st,
        N_MAPPED = nm,
        stringsAsFactors = FALSE
      )
    }
  }
  probe_df <- do.call(rbind, probe_rows)
  best_idx <- which.max(probe_df$N_MAPPED)
  id_from <- probe_df$KEYTYPE[best_idx]
  strategy <- probe_df$STRATEGY[best_idx]

  map_univ <- scq_try_map_strategy(feats, org.Hs.eg.db::org.Hs.eg.db, id_from, strategy = strategy)
  common_cols <- intersect(
    c("ENTREZID", "SYMBOL", "ALIAS", "ENSEMBL", "ENSEMBLTRANS", "ENSEMBLPROT"),
    colnames(map_univ)
  )
  map_univ <- unique(map_univ[, common_cols, drop = FALSE])
  universe_entrez <- unique(stats::na.omit(map_univ$ENTREZID))

  if (length(universe_entrez) < 0.05 * length(unique(feats))) {
    union_maps <- list(map_univ)
    for (kt in intersect(c("SYMBOL", "ALIAS", "ENSEMBL", "ENSEMBLTRANS", "ENSEMBLPROT"), kt_cand)) {
      for (st in c("raw", "title", "upper", "auto")) {
        tmp <- scq_try_map_strategy(feats, org.Hs.eg.db::org.Hs.eg.db, kt, strategy = st)
        if (is.data.frame(tmp) && nrow(tmp)) {
          tmp <- tmp[, intersect(
            c("ENTREZID", "SYMBOL", "ALIAS", "ENSEMBL", "ENSEMBLTRANS", "ENSEMBLPROT"),
            colnames(tmp)
          ), drop = FALSE]
          union_maps[[paste(kt, st, sep = ":")]] <- tmp
        }
      }
    }
    union_df <- do.call(rbind, union_maps)
    union_df <- unique(union_df[!is.na(union_df$ENTREZID) & nzchar(union_df$ENTREZID), , drop = FALSE])
    universe_entrez <- unique(union_df$ENTREZID)
  }

  genes_pool <- unique(top_per_cluster$gene)
  map_top <- scq_try_map_strategy(genes_pool, org.Hs.eg.db::org.Hs.eg.db, id_from, strategy = strategy)
  if (!nrow(map_top) || sum(!is.na(map_top$ENTREZID)) < 3) {
    union_top <- list(map_top)
    for (kt in intersect(c("SYMBOL", "ALIAS", "ENSEMBL", "ENSEMBLTRANS", "ENSEMBLPROT"), kt_cand)) {
      for (st in c("raw", "title", "upper", "auto")) {
        tmp <- scq_try_map_strategy(genes_pool, org.Hs.eg.db::org.Hs.eg.db, kt, strategy = st)
        if (is.data.frame(tmp) && nrow(tmp)) {
          union_top[[paste(kt, st, sep = ":")]] <- tmp
        }
      }
    }
    map_top <- unique(do.call(rbind, lapply(union_top, function(x) {
      x[, intersect(
        c("ENTREZID", "SYMBOL", "ALIAS", "ENSEMBL", "ENSEMBLTRANS", "ENSEMBLPROT"),
        colnames(x)
      ), drop = FALSE]
    })))
  }
  key_col <- intersect(c(id_from, "SYMBOL", "ALIAS", "ENSEMBL", "ENSEMBLTRANS", "ENSEMBLPROT"),
    colnames(map_top)
  )
  key_col <- if (length(key_col)) key_col[1] else "SYMBOL"
  key2ent <- setNames(map_top$ENTREZID, as.character(map_top[[key_col]]))

  go_results <- list()
  clusters <- heat$celltype_order
  min_needed <- min_genes_for_enrich
  if (min_needed > 5 && length(universe_entrez) < 2000) {
    min_needed <- max(5, floor(min_needed * 0.6))
  }

  for (ct in clusters) {
    g_keys <- unique(top_per_cluster$gene[top_per_cluster$cluster == ct])
    g_entrez <- unique(stats::na.omit(key2ent[as.character(g_keys)]))
    if (length(g_entrez) < min_needed || length(universe_entrez) < 100) {
      go_results[[ct]] <- data.frame(
        celltype = ct,
        ID = NA_character_,
        Description = "(no significant BP term)",
        GeneRatio = NA_character_,
        Count = 0L,
        p.adjust = 1,
        stringsAsFactors = FALSE
      )
      next
    }
    enr <- clusterProfiler::enrichGO(
      gene = g_entrez,
      OrgDb = org.Hs.eg.db::org.Hs.eg.db,
      keyType = "ENTREZID",
      ont = "BP",
      pAdjustMethod = "BH",
      universe = universe_entrez,
      qvalueCutoff = 0.2,
      pvalueCutoff = 0.05,
      readable = TRUE
    )
    df_full <- as.data.frame(enr)
    if (!nrow(df_full)) {
      go_results[[ct]] <- data.frame(
        celltype = ct,
        ID = NA_character_,
        Description = "(no significant BP term)",
        GeneRatio = NA_character_,
        Count = 0L,
        p.adjust = 1,
        stringsAsFactors = FALSE
      )
      next
    }
    n_sig <- sum(df_full$p.adjust <= 0.20, na.rm = TRUE)
    if (n_sig == 0) {
      df_full <- df_full[order(df_full$pvalue, -df_full$Count), , drop = FALSE]
      df_full <- utils::head(df_full, 1L)
    } else {
      df_full <- df_full[df_full$p.adjust <= 0.20, , drop = FALSE]
    }
    df_full$celltype <- ct
    df_full <- df_full[order(df_full$p.adjust, -df_full$Count), , drop = FALSE]
    go_results[[ct]] <- utils::head(
      df_full[, c("celltype", "ID", "Description", "GeneRatio", "p.adjust", "Count")],
      top_go_per_cluster
    )
  }

  go_cell <- do.call(rbind, go_results)
  go_cell <- as.data.frame(go_cell, stringsAsFactors = FALSE)
  go_cell$GeneRatio_decimal <- parse_ratio_decimal(go_cell$GeneRatio)
  go_cell$nl10FDR <- -log10(go_cell$p.adjust)
  go_cell$celltype <- factor(go_cell$celltype, levels = heat$celltype_order)
  go_cell$placeholder <- is.na(go_cell$ID) | grepl("^\\(no significant", go_cell$Description)
  go_cell <- go_cell[order(as.integer(go_cell$celltype), -go_cell$nl10FDR, -go_cell$Count), , drop = FALSE]
  go_cell$display_order <- seq_len(nrow(go_cell))
  go_cell$term_label <- paste0(as.character(go_cell$celltype), ": ", go_cell$Description)

  cl_vec <- as.character(obj$seurat_clusters)
  ct_vec <- as.character(obj@meta.data[[group.by]])
  maj_ct <- tapply(ct_vec, cl_vec, function(v) {
    v <- v[!is.na(v) & nzchar(v)]
    if (!length(v)) {
      NA_character_
    } else {
      names(sort(table(v), decreasing = TRUE))[1]
    }
  })

  id_backup <- Seurat::Idents(obj)
  on.exit({
    try(Seurat::Idents(obj) <- id_backup, silent = TRUE)
  }, add = TRUE)
  Seurat::Idents(obj) <- "seurat_clusters"
  cl_markers_all <- Seurat::FindAllMarkers(
    obj,
    assay = assay,
    only.pos = TRUE,
    logfc.threshold = 0.25,
    min.pct = 0.2,
    verbose = FALSE
  )

  cl_top <- lapply(split(cl_markers_all, cl_markers_all$cluster), function(df) {
    lfc_col <- intersect(c("avg_log2FC", "avg_logFC"), colnames(df))[1]
    if (is.na(lfc_col)) {
      lfc_col <- setdiff(colnames(df), c("gene", "cluster"))[1]
    }
    p_col <- if ("p_val_adj" %in% colnames(df)) {
      "p_val_adj"
    } else if ("p_val" %in% colnames(df)) {
      "p_val"
    } else {
      NULL
    }
    ord <- if (!is.null(p_col)) {
      order(-df[[lfc_col]], df[[p_col]], na.last = TRUE)
    } else {
      order(-df[[lfc_col]], na.last = TRUE)
    }
    head(df$gene[ord], top_n)
  })

  mat_u <- tryCatch(
    Seurat::GetAssayData(obj, assay = assay, layer = "data"),
    error = function(e) NULL
  )
  if (is.null(mat_u)) {
    mat_u <- Seurat::GetAssayData(obj, assay = assay, layer = "counts")
  }
  all_symbols <- rownames(mat_u)
  sym2ent_all <- suppressWarnings(
    clusterProfiler::bitr(
      all_symbols,
      fromType = "SYMBOL",
      toType = "ENTREZID",
      OrgDb = org.Hs.eg.db::org.Hs.eg.db
    )
  )
  universe_entrez2 <- unique(stats::na.omit(sym2ent_all$ENTREZID))

  genes_union <- unique(unlist(cl_top, use.names = FALSE))
  key2ent2 <- key2ent
  missing <- setdiff(genes_union, names(key2ent2))
  if (length(missing)) {
    add <- scq_try_map_strategy(missing, org.Hs.eg.db::org.Hs.eg.db, id_from, strategy = strategy)
    if (is.data.frame(add) && nrow(add) && "ENTREZID" %in% names(add)) {
      kcol <- intersect(c(id_from, "SYMBOL", "ALIAS", "ENSEMBL", "ENSEMBLTRANS", "ENSEMBLPROT"),
        colnames(add)
      )
      kcol <- if (length(kcol)) kcol[1] else "SYMBOL"
      add_map <- setNames(add$ENTREZID, as.character(add[[kcol]]))
      add_map <- add_map[setdiff(names(add_map), names(key2ent2))]
      key2ent2 <- c(key2ent2, add_map)
    }
  }

  all_clusters <- as.character(sort(unique(obj$seurat_clusters)))
  iter_clusters <- unique(c(all_clusters, names(cl_top)))
  go_rows <- list()
  for (cl in iter_clusters) {
    genes <- unique(cl_top[[cl]])
    if (is.null(genes) || !length(genes)) {
      go_rows[[cl]] <- data.frame(
        cluster = cl,
        ID = NA_character_,
        Description = "(no significant BP term)",
        GeneRatio = NA_character_,
        p.adjust = NA_real_,
        Count = 0L,
        stringsAsFactors = FALSE
      )
      next
    }
    g_entrez <- unique(stats::na.omit(unname(key2ent2[as.character(genes)])))
    if (!length(g_entrez)) {
      tmp <- scq_try_map_strategy(genes, org.Hs.eg.db::org.Hs.eg.db, id_from, strategy = strategy)
      if (is.data.frame(tmp) && nrow(tmp) && "ENTREZID" %in% names(tmp)) {
        g_entrez <- unique(stats::na.omit(tmp$ENTREZID))
      }
    }
    if (!length(g_entrez) || length(universe_entrez2) < 100) {
      go_rows[[cl]] <- data.frame(
        cluster = cl,
        ID = NA_character_,
        Description = "(no significant BP term)",
        GeneRatio = NA_character_,
        p.adjust = NA_real_,
        Count = 0L,
        stringsAsFactors = FALSE
      )
      next
    }
    ego <- suppressMessages(clusterProfiler::enrichGO(
      gene = g_entrez,
      universe = universe_entrez2,
      OrgDb = org.Hs.eg.db::org.Hs.eg.db,
      keyType = "ENTREZID",
      ont = "BP",
      pAdjustMethod = "BH",
      readable = TRUE
    ))
    df <- as.data.frame(ego)
    if (!nrow(df)) {
      go_rows[[cl]] <- data.frame(
        cluster = cl,
        ID = NA_character_,
        Description = "(no significant BP term)",
        GeneRatio = NA_character_,
        p.adjust = NA_real_,
        Count = 0L,
        stringsAsFactors = FALSE
      )
      next
    }
    df <- df[order(df$p.adjust, -df$Count), , drop = FALSE]
    df$cluster <- cl
    go_rows[[cl]] <- utils::head(
      df[, c("cluster", "ID", "Description", "GeneRatio", "p.adjust", "Count")],
      top_go_per_cluster
    )
  }

  go_cl <- do.call(rbind, go_rows)
  go_cl <- as.data.frame(go_cl, stringsAsFactors = FALSE)
  go_cl$GeneRatio_decimal <- parse_ratio_decimal(go_cl$GeneRatio)
  go_cl$nl10FDR <- -log10(go_cl$p.adjust)
  cl_order <- sort(unique(as.integer(as.character(go_cl$cluster))))
  go_cl$cluster <- factor(as.character(go_cl$cluster), levels = as.character(cl_order))
  lab_celltype <- maj_ct[as.character(go_cl$cluster)]
  lab_celltype[is.na(lab_celltype) | !nzchar(lab_celltype)] <- "Unknown"
  go_cl$majority_celltype <- as.character(lab_celltype)
  ord <- order(as.integer(as.character(go_cl$cluster)), -go_cl$nl10FDR, -go_cl$Count)
  go_cl <- go_cl[ord, , drop = FALSE]
  go_cl$display_order <- seq_len(nrow(go_cl))
  go_cl$term_label <- paste0(
    "C", as.character(go_cl$cluster), ": ",
    go_cl$majority_celltype, " _ ",
    go_cl$Description
  )

  list(
    celltype_go = go_cell,
    cluster_go = go_cl,
    keytype = id_from,
    strategy = strategy,
    universe = length(universe_entrez),
    universe_cluster = length(universe_entrez2)
  )
}

message("[1/7] Loading PBMC object and velocity layers")
obj <- readRDS(file.path(repo_dir, "benchmarking/qc_pbmc/pbmc_seu_processed.rds"))
sp <- readRDS(file.path(repo_dir, "benchmarking/pbmc10k/testloom.rds"))

message("[2/7] Recreating report auto-annotation")
sig <- panglao_signatures(
  tsv = file.path(repo_dir, "inst/extdata/PanglaoDB_markers_27_Mar_2020.tsv"),
  species = "human",
  canonical_only = TRUE,
  ui_max = 0.20,
  tissue = "Immune system",
  min_genes = 5
)
obj <- auto_annotate(
  obj,
  species = "human",
  assay = "RNA",
  prefer = "marker_score",
  level = "cluster",
  marker_score_level = "cluster",
  custom_signatures = sig,
  marker_method = "module_score",
  unknown_min_genes = 2,
  unknown_min_margin = 0.01,
  unknown_min_top = 0.01
)

message("[3/7] Recomputing QC metrics used by the feature maps")
metrics <- calcQCmetrics(
  obj,
  assay = "RNA",
  slot = "counts",
  species = "human",
  spliced_obj = sp,
  spliced_assay = "spliced",
  spliced_layer = "counts",
  unspliced_obj = sp,
  unspliced_assay = "unspliced",
  unspliced_layer = "counts",
  add_to_meta = TRUE
)
obj <- attr(metrics, "obj")
if (!("is_doublet" %in% colnames(obj@meta.data))) {
  obj@meta.data$is_doublet <- FALSE
}
if (!("dbl_score" %in% colnames(obj@meta.data))) {
  obj@meta.data$dbl_score <- NA_real_
}
if (!("dbl_method" %in% colnames(obj@meta.data))) {
  obj@meta.data$dbl_method <- NA_character_
}

umap <- get_umap(obj)
md <- obj@meta.data[umap$cell, , drop = FALSE]

message("[4/7] Writing Figure 6A/B/F source tables")
fig6a <- cbind(
  umap,
  seurat_cluster = as.character(md$seurat_clusters),
  auto_celltype = as.character(md$auto_celltype),
  auto_celltype_source = as.character(md$auto_celltype_src),
  auto_celltype_confidence = as.numeric(md$auto_celltype_conf),
  stringsAsFactors = FALSE
)
write_csv_base(fig6a, file.path(out_dir, "Figure 6A.csv"))

metric_specs <- data.frame(
  metric = c(
    "nCount_QC", "nFeature_QC", "pctMT_QC", "MALAT1_frac_QC",
    "stress_score_QC", "intronic_frac_QC", "spliced_frac_QC", "unspliced_frac_QC"
  ),
  feature_order = seq_len(8),
  plot_title = c(
    "Total counts (log)", "Detected genes (log)", "% mito", "MALAT1 fraction",
    "Stress score", "Intronic fraction", "Spliced fraction", "Unspliced fraction"
  ),
  log10_transform = c(TRUE, TRUE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE),
  stringsAsFactors = FALSE
)

fig6b_parts <- lapply(seq_len(nrow(metric_specs)), function(i) {
  spec <- metric_specs[i, ]
  raw <- as.numeric(md[[spec$metric]])
  plot_value <- raw
  if (isTRUE(spec$log10_transform)) {
    plot_value <- log10(plot_value + 1)
  }
  plot_value <- winsorize(plot_value, c(0.02, 0.98))
  data.frame(
    cell = umap$cell,
    UMAP_1 = umap$UMAP_1,
    UMAP_2 = umap$UMAP_2,
    metric = spec$metric,
    feature_order = spec$feature_order,
    plot_title = spec$plot_title,
    raw_value = raw,
    plot_value = plot_value,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
})
fig6b <- do.call(rbind, fig6b_parts)
write_csv_base(fig6b, file.path(out_dir, "Figure 6B.csv"))

fig6f <- cbind(
  umap,
  seurat_cluster = as.character(md$seurat_clusters),
  auto_celltype = as.character(md$auto_celltype),
  qc_status = factor(as.character(md$qc_status), levels = c("keep", "borderline", "remove")),
  qc_score = if ("qc_score" %in% colnames(md)) as.numeric(md$qc_score) else NA_real_,
  retained_for_report = as.character(md$qc_status) %in% c("keep", "borderline"),
  stringsAsFactors = FALSE
)
write_csv_base(fig6f, file.path(out_dir, "Figure 6F.csv"))

message("[5/7] Recomputing heatmap and enrichment source values")
heat <- make_heatmap_source(obj, assay = "RNA", group.by = "auto_celltype", slot.use = "data")
write_csv_base(heat$heatmap_source, file.path(out_dir, "Figure 6C.csv"))

go <- go_sources(
  obj,
  heat,
  assay = "RNA",
  group.by = "auto_celltype",
  top_n = 100,
  min_genes_for_enrich = 10,
  top_go_per_cluster = 2
)

fig6d <- go$celltype_go[, c(
  "celltype", "display_order", "ID", "Description",
  "GeneRatio_decimal", "p.adjust", "nl10FDR"
), drop = FALSE]
names(fig6d)[names(fig6d) == "GeneRatio_decimal"] <- "GeneRatio"
write_csv_base(fig6d, file.path(out_dir, "Figure 6D.csv"))

fig6e <- go$cluster_go[, c(
  "cluster", "majority_celltype", "display_order", "ID", "Description",
  "GeneRatio_decimal", "p.adjust", "nl10FDR"
), drop = FALSE]
names(fig6e)[names(fig6e) == "GeneRatio_decimal"] <- "GeneRatio"
write_csv_base(fig6e, file.path(out_dir, "Figure 6E.csv"))

message("[6/7] Writing verification plots from source tables")
p6a <- ggplot(fig6a, aes(UMAP_1, UMAP_2, color = auto_celltype)) +
  geom_point(size = 0.32, alpha = 0.68) +
  coord_equal() +
  theme_classic(base_size = 16) +
  labs(title = "UMAP (auto_celltype)", color = "auto_celltype") +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1)))
ggsave(file.path(verify_dir, "Figure 6A_verify.png"), p6a, width = 9, height = 7, dpi = 200, bg = "white")

b_plots <- lapply(split(fig6b, fig6b$feature_order), function(d) {
  feature_plot_from_source(d, unique(d$plot_title))
})
p6b <- patchwork::wrap_plots(b_plots, ncol = 4, byrow = TRUE, guides = "keep")
ggsave(file.path(verify_dir, "Figure 6B_verify.png"), p6b, width = 20, height = 9.5, dpi = 180, bg = "white")

heat_plot_df <- heat$heatmap_source
heat_plot_df$gene <- factor(heat_plot_df$gene, levels = rev(rownames(heat$heatmap_matrix)))
heat_plot_df$auto_celltype <- factor(heat_plot_df$auto_celltype, levels = colnames(heat$heatmap_matrix))
p6c <- ggplot(heat_plot_df, aes(auto_celltype, gene, fill = average_expression)) +
  geom_tile() +
  scale_fill_viridis_c(name = "expr") +
  theme_minimal(base_size = 9) +
  labs(x = NULL, y = NULL) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    axis.text.y = element_text(size = 7)
  )
ggsave(file.path(verify_dir, "Figure 6C_verify.png"), p6c, width = 7.5, height = 11, dpi = 200, bg = "white")

plot_d <- go$celltype_go
plot_d$term_label <- factor(plot_d$term_label, levels = rev(plot_d$term_label))
p6d <- ggplot(plot_d, aes(nl10FDR, term_label, fill = celltype)) +
  geom_col() +
  theme_minimal(base_size = 10) +
  labs(x = expression(-log[10]("adj. P")), y = NULL, title = "Top GO BP term per cell type") +
  theme(legend.position = "none", panel.grid = element_blank())
ggsave(file.path(verify_dir, "Figure 6D_verify.png"), p6d, width = 8.5, height = nrow(plot_d) * 0.15 + 1, dpi = 200, bg = "white")

plot_e <- go$cluster_go
plot_e$term_label <- factor(plot_e$term_label, levels = rev(plot_e$term_label))
p6e <- ggplot(plot_e, aes(nl10FDR, term_label, fill = cluster)) +
  geom_col() +
  theme_minimal(base_size = 10) +
  labs(x = expression(-log[10]("adj. P")), y = NULL, title = "Top GO BP term per Seurat cluster") +
  theme(legend.position = "none", panel.grid = element_blank())
ggsave(file.path(verify_dir, "Figure 6E_verify.png"), p6e, width = 8.5, height = nrow(plot_e) * 0.15 + 1, dpi = 200, bg = "white")

p6f <- ggplot(fig6f, aes(UMAP_1, UMAP_2, color = qc_status)) +
  geom_point(size = 0.32, alpha = 0.68) +
  coord_equal() +
  theme_classic(base_size = 16) +
  scale_color_manual(values = c(keep = "#49A879", borderline = "#5DA5DA", remove = "#F05B61"), drop = FALSE) +
  labs(title = "UMAP (colored by QC status)", color = "qc_status") +
  guides(color = guide_legend(override.aes = list(size = 3, alpha = 1)))
ggsave(file.path(verify_dir, "Figure 6F_verify.png"), p6f, width = 8, height = 7, dpi = 200, bg = "white")

message("[7/7] Writing report component plots and verification summary")
qc_diagnostic_panel(
  obj,
  metrics,
  status_df = data.frame(qc_status = md$qc_status, row.names = rownames(md)),
  assay = "RNA",
  save_dir = component_dir,
  prefix = "qc",
  group_col = "auto_celltype"
)
qc_featuremaps(
  obj,
  metrics = metric_specs$metric,
  save_dir = component_dir,
  prefix = "qc",
  ncol = 4
)

summary_lines <- c(
  "Figure 6 report source-data regeneration",
  sprintf("Generated: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "",
  "Inputs:",
  file.path(repo_dir, "benchmarking/qc_pbmc/pbmc_seu_processed.rds"),
  file.path(repo_dir, "benchmarking/pbmc10k/testloom.rds"),
  "",
  "Annotation settings:",
  "species = human; tissue = Immune system; marker_score_level = cluster; marker_method = module_score",
  "",
  "Panel source files:",
  file.path(out_dir, "Figure 6A.csv"),
  file.path(out_dir, "Figure 6B.csv"),
  file.path(out_dir, "Figure 6C.csv"),
  file.path(out_dir, "Figure 6D.csv"),
  file.path(out_dir, "Figure 6E.csv"),
  file.path(out_dir, "Figure 6F.csv"),
  "",
  "Verification plots:",
  file.path(verify_dir, "Figure 6A_verify.png"),
  file.path(verify_dir, "Figure 6B_verify.png"),
  file.path(verify_dir, "Figure 6C_verify.png"),
  file.path(verify_dir, "Figure 6D_verify.png"),
  file.path(verify_dir, "Figure 6E_verify.png"),
  file.path(verify_dir, "Figure 6F_verify.png"),
  "",
  sprintf("Figure 6A cells: %d", nrow(fig6a)),
  sprintf("Figure 6B source rows: %d (%d metrics x %d cells)", nrow(fig6b), nrow(metric_specs), nrow(fig6a)),
  sprintf("Figure 6C matrix: %d genes x %d cell types", nrow(heat$heatmap_matrix), ncol(heat$heatmap_matrix)),
  sprintf("Figure 6D enrichment rows: %d", nrow(fig6d)),
  sprintf("Figure 6E enrichment rows: %d", nrow(fig6e)),
  sprintf("Figure 6F QC counts: %s", paste(names(table(fig6f$qc_status)), as.integer(table(fig6f$qc_status)), sep = "=", collapse = ", ")),
  sprintf("GO mapping keytype: %s; strategy: %s; universe=%d; cluster_universe=%d",
    go$keytype, go$strategy, go$universe, go$universe_cluster
  )
)
writeLines(summary_lines, file.path(out_dir, "Figure_6_README.txt"), useBytes = TRUE)

message(paste(summary_lines, collapse = "\n"))
