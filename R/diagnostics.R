#' Ensure PCA/UMAP are present for nice plots (Seurat v5)
#' @keywords internal
.ensure_visual_embeddings <- function(seu, assay = "RNA", pc_dims = 1:20,
                                      debug = getOption("scQCenrich.debug", FALSE)) {
  stopifnot(inherits(seu, "Seurat"))
  Seurat::DefaultAssay(seu) <- assay

  pca_ok  <- "pca"  %in% names(seu@reductions)
  umap_ok <- "umap" %in% names(seu@reductions)

  if (!pca_ok || !umap_ok) {
    if (debug) message("[scQCenrich-debug] NormalizeData(): creating non-empty 'data' layer for diagnostics.")
    seu <- Seurat::NormalizeData(seu, verbose = FALSE)
    seu <- Seurat::FindVariableFeatures(seu, nfeatures = 2000, verbose = FALSE)
    hvgs <- Seurat::VariableFeatures(seu)

    if (debug) message("[scQCenrich-debug] FindVariableFeatures/ScaleData/RunPCA(): creating PCA.")
    if (length(hvgs)) {
      seu <- Seurat::ScaleData(seu, features = hvgs, verbose = FALSE)
      seu <- Seurat::RunPCA(seu, features = hvgs, npcs = max(pc_dims), verbose = FALSE)
    } else {
      seu <- Seurat::ScaleData(seu, verbose = FALSE)
      seu <- Seurat::RunPCA(seu, npcs = max(pc_dims), verbose = FALSE)
    }

    if (debug) message("[scQCenrich-debug] RunUMAP(): creating UMAP with ", max(pc_dims), " PCs.")
    nd <- min(max(pc_dims), ncol(Seurat::Embeddings(seu, "pca")))
    seu <- Seurat::RunUMAP(seu, dims = 1:nd, verbose = FALSE, seed.use = .umap_seed())
  }
  seu
}





# at top-level helpers, pick a default seed (can be overridden by options())
.umap_seed <- function() getOption("scQCenrich.umap_seed", 1337L)


