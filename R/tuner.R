#' Grid-search to preserve small clusters (k.param & resolution)
#' @param obj Seurat object
#' @param dims_use Integer vector of PC indices
#' @param k_vals Integer vector of k.param candidates
#' @param res_vals Numeric vector of resolution candidates
#' @param algorithm Integer; clustering algorithm (e.g., Leiden=4)
#' @param prune Numeric; prune.SNN
#' @param small_frac_target Length-2 numeric; desired min cluster fraction range
#' @param seed Integer random seed
#' @return List with components $stats (grid) and $pick (chosen row)
#' @export

tune_clustering <- function(
    obj,
    dims_use = 1:30,
    k_vals   = c(10, 12, 15, 20),
    res_vals = seq(0.2, 1.6, by = 0.2),
    algorithm = 4,     # Leiden
    prune = 0.1,
    small_frac_target = c(0.005, 0.05),  # 0.5%-5%
    seed = 1
){
  stopifnot(inherits(obj, "Seurat"))
  set.seed(seed)
  n <- ncol(obj)

  if (is.null(obj[["pca"]])) {
    obj <- Seurat::NormalizeData(obj, verbose = FALSE)
    obj <- Seurat::FindVariableFeatures(obj, nfeatures = 2000, verbose = FALSE)
    obj <- Seurat::ScaleData(obj, verbose = FALSE)
    obj <- Seurat::RunPCA(obj, npcs = max(dims_use), verbose = FALSE)
  }

  results <- list()
  for (k in k_vals) {
    gname <- paste0("snn_k", k)
    obj <- Seurat::FindNeighbors(obj, dims = dims_use, k.param = k, prune.SNN = prune,
                                 graph.name = gname, verbose = FALSE)
    for (res in res_vals) {
      obj <- Seurat::FindClusters(obj, graph.name = gname, resolution = res,
                                  algorithm = algorithm, random.seed = seed, verbose = FALSE)
      cl  <- obj$seurat_clusters
      tab <- sort(table(cl), decreasing = TRUE)
      min_sz   <- min(tab)
      min_frac <- as.numeric(min_sz) / n

      mean_sil <- NA_real_
      if (requireNamespace("cluster", quietly = TRUE)) {
        emb <- Seurat::Embeddings(obj, "pca")[, dims_use, drop = FALSE]
        sil <- try(cluster::silhouette(as.integer(cl), stats::dist(emb)), silent = TRUE)
        if (!inherits(sil, "try-error")) mean_sil <- mean(sil[,3])
      }

      results[[length(results) + 1]] <- data.frame(
        k = k, res = res, ncl = length(tab),
        min_size = as.numeric(min_sz), min_frac = min_frac,
        mean_sil = mean_sil, stringsAsFactors = FALSE
      )
    }
  }
  stats <- do.call(rbind, results)

  lo <- small_frac_target[1]; hi <- small_frac_target[2]
  cand <- subset(stats, min_frac >= lo & min_frac <= hi)
  if (!nrow(cand)) cand <- stats[order(abs(stats$min_frac - mean(c(lo,hi))), stats$res), ][1, , drop = FALSE]
  cand <- cand[order(cand$res, -ifelse(is.na(cand$mean_sil), -Inf, cand$mean_sil), cand$k), ][1, , drop = FALSE]

  list(stats = stats, pick = cand)
}
#' Apply chosen k/res to an object (returns updated object)
#'
#' @param obj Seurat object to update
#' @param dims_use Integer vector of PC indices used to build the graph
#' @param k_best Integer; chosen k.param
#' @param res_best Numeric; chosen resolution
#' @param prune Numeric; prune.SNN passed to FindNeighbors
#' @param algorithm Integer; clustering algorithm (e.g., 4 = Leiden)
#' @param seed Integer random seed
#' @return Updated Seurat object
#' @export

apply_tuned_clustering <- function(obj, dims_use, k_best, res_best, prune = 0.1, algorithm = 4, seed = 1){
  gname <- paste0("snn_k", k_best)
  obj <- Seurat::FindNeighbors(obj, dims = dims_use, k.param = k_best, prune.SNN = prune,
                               graph.name = gname, verbose = FALSE)
  obj <- Seurat::FindClusters(obj, graph.name = gname, algorithm = algorithm,
                              resolution = res_best, random.seed = seed, verbose = FALSE)
  obj
}
