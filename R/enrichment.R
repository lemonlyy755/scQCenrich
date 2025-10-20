#' Prepare ranked gene stats for dropped vs kept (logCPM mean diff)
#' @param obj Seurat object
#' @param status_df Data frame with rownames=cells and 'qc_status'
#' @param assay Assay name to use
#' @return Named numeric vector of statistics (higher = up in dropped)
#' @export

prep_flagged_ranks <- function(obj, status_df, assay = "RNA"){
  mat <- Seurat::GetAssayData(obj, layer  = "counts", assay = assay)
  if (ncol(mat) == 0L) stop("prep_flagged_ranks(): counts matrix has 0 columns.")

  if (!inherits(mat, "dgCMatrix")) mat <- methods::as(mat, "dgCMatrix")
  keep <- intersect(colnames(mat), rownames(status_df))
  if (!length(keep)) stop("No overlap between object and status.")
  logcpm <- .log_cpm(mat[,keep,drop=FALSE])
  flag <- status_df[keep, "qc_status"] == "remove"

  if (all(flag) || !any(flag)) stop("Need both drop and keep cells.")
  mu_drop <- Matrix::rowMeans(logcpm[,flag,drop=FALSE])
  mu_keep <- Matrix::rowMeans(logcpm[,!flag,drop=FALSE])
  stat <- (mu_drop - mu_keep); names(stat) <- rownames(logcpm)
  sort(stat[is.finite(stat)], decreasing = TRUE)
}

#' Hallmark enrichment (ORA + GSEA) for dropped cells
#' @param ranks Named numeric vector of gene statistics
#' @param species 'mouse' or 'human'
#' @param top_n_up Number of top up-genes to use for ORA
#' @return List with data.frames: $ora and $gsea (or NULLs if unavailable)
#' @export

run_enrichment <- function(ranks, species = c("mouse","human"), top_n_up = 200){
  species <- match.arg(species)
  if (!requireNamespace("msigdbr", quietly = TRUE) ||
      !requireNamespace("clusterProfiler", quietly = TRUE) ||
      !requireNamespace("fgsea", quietly = TRUE)) {
    warning("Install msigdbr, clusterProfiler, fgsea for enrichment.")
    return(list(ora=NULL, gsea=NULL))
  }

  sp <- if (species == "mouse") "Mus musculus" else "Homo sapiens"

  # msigdbr v10+ uses `collection` (not `category`)
  msig <- msigdbr::msigdbr(species = sp, collection = "H")
  TERM2GENE <- msig[, c("gs_name","gene_symbol")]

  up_genes <- names(ranks)[ranks > 0][1:min(top_n_up, sum(ranks > 0))]
  ora <- try(clusterProfiler::enricher(up_genes, TERM2GENE = TERM2GENE), silent = TRUE)

  hallmark <- split(msig$gene_symbol, f = gsub("^HALLMARK_", "", msig$gs_name))
  gsea <- try(fgsea::fgseaMultilevel(pathways = hallmark, stats = ranks,
                                     minSize = 10, maxSize = 1000), silent = TRUE)

  gsea_df <- NULL
  if (!inherits(gsea, "try-error")) {
    gsea_df <- as.data.frame(gsea)

    list_cols <- vapply(gsea_df, is.list, logical(1))
    if (any(list_cols)) {
      for (nm in names(gsea_df)[list_cols]) {
        gsea_df[[nm]] <- vapply(gsea_df[[nm]], function(x) paste(x, collapse = ";"), "")
      }
    }
    gsea_df <- gsea_df[order(gsea_df$padj, -abs(gsea_df$NES)), ]
  }

  list(
    ora  = if (!inherits(ora,  "try-error")) as.data.frame(ora)  else NULL,
    gsea = gsea_df
  )
}

#' Removed cells clustering + enrichment (used by report)
#' @param obj Seurat object
#' @param status_df Data frame with 'qc_status'
#' @param save_dir Output directory for PNG/CSV
#' @param species 'mouse' or 'human'
#' @param assay Assay name
#' @return Invisibly returns NULL; writes files to disk
#' @export

removed_cells_analysis <- function(obj, status_df, save_dir="qc_outputs", species=c("mouse","human"), assay="RNA"){
  species <- match.arg(species)
  drops <- rownames(status_df)[status_df$qc_status == "remove"]

  if (length(drops) < 50) return(invisible(NULL))
  seu <- subset(obj, cells = drops)
  Seurat::DefaultAssay(seu) <- assay
  seu <- Seurat::NormalizeData(seu, verbose=FALSE)
  seu <- Seurat::FindVariableFeatures(seu, selection.method = "vst", nfeatures = 2000, verbose=FALSE)
  seu <- Seurat::ScaleData(seu, verbose=FALSE)
  seu <- Seurat::RunPCA(seu, npcs = 30, verbose=FALSE)
  seu <- Seurat::RunUMAP(seu, dims = 1:20, verbose=FALSE)
  seu <- Seurat::FindNeighbors(seu, dims = 1:20, verbose=FALSE)
  seu <- Seurat::FindClusters(seu, resolution = 0.4, verbose=FALSE)
  gr <- try(Seurat::DimPlot(seu, reduction="umap", group.by="seurat_clusters", label=TRUE) +
              ggplot2::ggtitle("Removed cells (UMAP, Seurat clusters)"),
            silent = TRUE)
  if (!inherits(gr, "try-error")) .gg_save_safe(gr, file.path(save_dir, "removed_cells_umap.png"), 8, 6)


  if (!requireNamespace("msigdbr", quietly = TRUE) ||
      !requireNamespace("clusterProfiler", quietly = TRUE)) return(invisible(NULL))

  Seurat::Idents(seu) <- "seurat_clusters"
  markers <- try(Seurat::FindAllMarkers(seu, only.pos=TRUE, logfc.threshold=0.25, min.pct=0.1, verbose=FALSE), silent=TRUE)
  if (inherits(markers, "try-error")) return(invisible(NULL))
  sp <- if (species=="mouse") "Mus musculus" else "Homo sapiens"
  msig <- msigdbr::msigdbr(species = sp, collection = "H")
  TERM2GENE <- msig[, c("gs_name","gene_symbol")]
  enr_all <- list()
  for (cl in sort(unique(markers$cluster))) {
    genes <- head(markers$gene[markers$cluster==cl], 200)
    erf <- try(clusterProfiler::enricher(genes, TERM2GENE = TERM2GENE), silent=TRUE)
    if (!inherits(erf, "try-error")) enr_all[[as.character(cl)]] <- as.data.frame(erf)
  }
  if (length(enr_all)) {
    all_tab <- data.frame(
      cluster = rep(names(enr_all), sapply(enr_all, nrow)),
      do.call(rbind, lapply(enr_all, function(x) x[, c("ID","Description","p.adjust","GeneRatio","Count"), drop=FALSE])),
      check.names = FALSE
    )
    utils::write.csv(all_tab, file.path(save_dir, "removed_cells_enrichment.csv"), row.names = FALSE)
  }
  invisible(NULL)
}
