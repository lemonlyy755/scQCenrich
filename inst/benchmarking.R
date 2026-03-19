


# install if missing ------------------------------------------------------

# 1. Increase timeout for large Bioconductor/GitHub downloads
options(timeout = 600)

# 2. Define the packages by source
cran_pkgs <- c("Seurat", "Matrix", "dplyr", "tidyr", "ggplot2", "cluster", "lme4", "remotes")
bioc_pkgs <- c("scater", "SingleCellExperiment", "scuttle", "miQC", "scDblFinder",
               "celda", "SingleR", "celldex")
github_pkgs <- c("DropletQC" = "powellgenomicslab/DropletQC")

# 3. Ensure BiocManager and remotes are installed first
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")

# 4. Install missing CRAN packages
missing_cran <- cran_pkgs[!(cran_pkgs %in% installed.packages()[, "Package"])]
if (length(missing_cran)) install.packages(missing_cran)

# 5. Install missing Bioconductor packages
missing_bioc <- bioc_pkgs[!(bioc_pkgs %in% installed.packages()[, "Package"])]
if (length(missing_bioc)) BiocManager::install(missing_bioc, ask = FALSE)

# 6. Install missing GitHub packages (specifically DropletQC)
if (!requireNamespace("DropletQC", quietly = TRUE)) {
  remotes::install_github(github_pkgs["DropletQC"])
}

# 7. Load all libraries
all_pkgs <- c(cran_pkgs, bioc_pkgs, names(github_pkgs))
lapply(all_pkgs, library, character.only = TRUE)


# load lib ----------------------------------------------------------------



library(scater) # addPerCellQCMetrics
library(Seurat)
library(SingleCellExperiment)
library(scuttle)
library(miQC) # mixtureModel, filterCells
library(scDblFinder) # doublets
library(celda) # DecontX
library(SingleR)
library(celldex) # biology annotation
library(Matrix)
library(dplyr)
library(tidyr)
library(ggplot2)
library(cluster) # silhouette
library(lme4) # mixed model
library(DropletQC) # dropletqc

options(timeout = 600)


# Data path-----

raw_counts_dir <- "path/to/filtered_feature_bc_matrix/"
spliced_rds <- "path/to/testloom.rds"
# human or mouse
organism <- "human"
organ <- c("Lung", "Blood")


# preprocessing -------------------------------------------------------
seurat_std <- CreateSeuratObject(counts = Read10X(raw_counts_dir), project = "QC_compare")

seurat_sp <- readRDS(spliced_rds)

# Helper: normalize cell barcodes to a comparable format
.normalize_cells <- function(cells) {
  # drop sample prefix if present
  cells <- sub(".*:", "", cells)
  # drop trailing x if present
  cells <- sub("x$", "", cells)
  # drop trailing -1 if present
  cells <- sub("-1$", "", cells)
  cells
}

# Suppose you have:
#   seurat_std  = your "standard" Seurat object (RNA assay)
#   seurat_sp   = your velocity-style Seurat object (spliced/unspliced)

# Normalize cell IDs
std_ids <- .normalize_cells(colnames(seurat_std))
sp_ids <- .normalize_cells(colnames(seurat_sp))

# Find intersection
common <- intersect(std_ids, sp_ids)
message("Found ", length(common), " overlapping cells.")

# Sample cells if needed (or fewer if less available)
set.seed(123) # for reproducibility
keep_ids <- sample(common, min(50000, length(common)))

# Map back to original names in each object
std_keep <- colnames(seurat_std)[.normalize_cells(colnames(seurat_std)) %in% keep_ids]
sp_keep <- colnames(seurat_sp)[.normalize_cells(colnames(seurat_sp)) %in% keep_ids]

# Subset
seurat_std_small <- subset(seurat_std, cells = std_keep)
seurat_sp_small <- subset(seurat_sp, cells = sp_keep)

# Check dimensions
dim(seurat_std_small) # genes x 500 (or fewer)
dim(seurat_sp_small)

seurat_std_small <- NormalizeData(seurat_std_small)
seurat_std_small <- FindVariableFeatures(seurat_std_small)
seurat_std_small <- ScaleData(seurat_std_small)
seurat_std_small <- RunPCA(seurat_std_small)
seurat_std_small <- RunUMAP(seurat_std_small, dims = 1:20, verbose = FALSE, seed.use = 1337)
seurat_std_small <- FindNeighbors(seurat_std_small)
seurat_std_small <- FindClusters(seurat_std_small, resolution = 0.4)


# metric functions --------------------------------------------------------
# ===========================
# Helpers (species-agnostic mito detection)
# ===========================
.mito_patterns_default <- c("^MT-", "^mt-", "^Mt-") # human, mouse, mixed
.mito_regex <- function(patterns) paste0("(", paste(patterns, collapse = "|"), ")")

.get_mito_features <- function(seu, mito_patterns = NULL) {
  if (is.null(mito_patterns)) mito_patterns <- .mito_patterns_default
  genes <- rownames(seu)
  if (is.null(genes)) {
    return(character(0))
  }
  genes[grepl(.mito_regex(mito_patterns), genes)]
}

.get_mito_index_sce <- function(sce, mito_patterns = NULL) {
  if (is.null(mito_patterns)) mito_patterns <- .mito_patterns_default
  rn <- rownames(sce)
  if (is.null(rn)) {
    return(rep(FALSE, nrow(sce)))
  }
  grepl(.mito_regex(mito_patterns), rn)
}

# ===========================
# 1) Utility: safe add meta by barcode  (unchanged)
# ===========================
add_meta <- function(seu, df, suffix = NULL) {
  stopifnot(all(rownames(df) %in% colnames(seu)))
  if (!is.null(suffix)) colnames(df) <- paste0(colnames(df), "_", suffix)
  df <- df[colnames(seu), , drop = FALSE]
  Seurat::AddMetaData(seu, df)
}

# ---- META-ONLY MERGE (donor -> target)  (minor namespace tighten; no behavior change) ----
merge_meta_only <- function(target, donor, cols = NULL, suffix = "_donor",
                            normalizer = .normalize_cells, verbose = TRUE) {
  stopifnot(inherits(target, "Seurat"), inherits(donor, "Seurat"))

  # maps of original -> normalized IDs
  map <- function(obj) {
    data.frame(
      orig = colnames(obj),
      norm = normalizer(colnames(obj)),
      stringsAsFactors = FALSE
    )
  }
  mT <- map(target)
  mD <- map(donor)

  # if normalization creates duplicate IDs, fall back to exact IDs
  if (any(duplicated(mT$norm)) || any(duplicated(mD$norm))) {
    if (verbose) message("Duplicates after normalization; falling back to exact barcode match.")
    mT$norm <- mT$orig
    mD$norm <- mD$orig
  }

  # intersection and index vectors
  common <- intersect(mT$norm, mD$norm)
  if (!length(common)) stop("No overlapping cells between target and donor after alignment.")
  idxT <- match(common, mT$norm)
  idxD <- match(common, mD$norm)

  # donor metadata (optionally subset columns)
  donor_meta <- donor@meta.data
  if (!is.null(cols)) {
    missing_cols <- setdiff(cols, colnames(donor_meta))
    if (length(missing_cols)) warning("Skipping missing donor columns: ", paste(missing_cols, collapse = ", "))
    cols <- intersect(cols, colnames(donor_meta))
    donor_meta <- donor_meta[, cols, drop = FALSE]
  }

  # align rows to target’s order; fill NA where target has non-overlapping cells
  aligned <- donor_meta[mD$orig[idxD], , drop = FALSE]
  rownames(aligned) <- mT$orig[idxT]
  out <- matrix(NA,
                nrow = ncol(target), ncol = ncol(aligned),
                dimnames = list(colnames(target), colnames(aligned))
  )
  out[rownames(aligned), ] <- as.matrix(aligned)
  out <- as.data.frame(out, check.names = FALSE)

  # handle name collisions
  collide <- intersect(colnames(target@meta.data), colnames(out))
  if (length(collide)) colnames(out)[match(collide, colnames(out))] <- paste0(collide, suffix)

  # add to target meta (by barcode rownames)
  target <- Seurat::AddMetaData(target, out)
  if (verbose) {
    message("Merged ", ncol(out), " column(s) onto target for ", length(common), " overlapping cells.")
  }
  return(target)
}

# ===========================
# 2) Baseline: Seurat thresholds (species-agnostic via features list)
# ===========================
baseline_seurat <- function(seu, min_genes = 200, max_genes = 10000, max_mito = 10, mito_patterns = NULL) {
  mito_features <- .get_mito_features(seu, mito_patterns)
  if (length(mito_features)) {
    seu[["percent.mt"]] <- Seurat::PercentageFeatureSet(seu, features = mito_features)
  } else {
    warning("No mitochondrial features matched; setting percent.mt to NA.")
    seu[["percent.mt"]] <- NA_real_
  }
  keep <- with(
    seu@meta.data,
    nFeature_RNA > min_genes & nFeature_RNA < max_genes &
      (is.na(percent.mt) | percent.mt < max_mito)
  )
  add_meta(seu, data.frame(keep_seurat = keep, row.names = colnames(seu)))
}

# ===========================
# 3) Baseline: scater/scuttle outlier rules (species-agnostic)
# ===========================
baseline_scater <- function(seu, mito_patterns = NULL) {
  sce <- as.SingleCellExperiment(seu)
  is_mito <- .get_mito_index_sce(sce, mito_patterns)
  sce <- scuttle::addPerCellQCMetrics(sce, subsets = list(mito = is_mito))
  low_lib <- scuttle::isOutlier(SummarizedExperiment::colData(sce)$sum, nmads = 3, type = "lower", log = TRUE)
  low_feat <- scuttle::isOutlier(SummarizedExperiment::colData(sce)$detected, nmads = 3, type = "lower", log = TRUE)
  high_mito <- scuttle::isOutlier(SummarizedExperiment::colData(sce)$subsets_mito_percent, nmads = 3, type = "higher")
  keep <- !(low_lib | low_feat | high_mito)
  add_meta(seu, data.frame(
    keep_scater = keep,
    subsets_mito_percent = SummarizedExperiment::colData(sce)$subsets_mito_percent,
    row.names = colnames(seu)
  ))
}

# ===========================
# 4) Baseline: miQC (robust + species-agnostic; defaults unchanged)
# ===========================
baseline_miQC <- function(seu, mito_patterns = NULL, debug = TRUE) {
  stopifnot(inherits(seu, "Seurat"))
  sce <- as.SingleCellExperiment(seu)
  # ensure counts exist
  if (is.null(SummarizedExperiment::assay(sce, "counts"))) {
    SummarizedExperiment::assay(sce, "counts") <- Seurat::GetAssayData(seu, layer = "counts")
  }

  # Prepare required metrics (detected & subsets_mito_percent) in a species-agnostic way
  is_mito <- .get_mito_index_sce(sce, mito_patterns)
  sce <- scater::addPerCellQC(sce, subsets = list(mito = rownames(sce)[is_mito])) # adds 'detected' & 'subsets_mito_percent'

  det <- as.numeric(SummarizedExperiment::colData(sce)$detected)
  mito_pct <- as.numeric(SummarizedExperiment::colData(sce)$subsets_mito_percent) # 0–100

  # Guard against non-finite/degenerate inputs to flexmix
  bad <- !is.finite(det) | !is.finite(mito_pct)
  if (any(bad)) {
    if (isTRUE(debug)) message(sprintf("[miQC] dropping %d cells with non-finite QC metrics before fit", sum(bad)))
    sce <- sce[, !bad, drop = FALSE]
  }
  mito_now <- as.numeric(SummarizedExperiment::colData(sce)$subsets_mito_percent)
  if (length(unique(mito_now)) <= 1L) {
    # Avoid zero-variance crashes; if still degenerate, skip miQC gracefully
    if (isTRUE(debug)) message("[miQC] degenerate mitochondrial metric; skipping model and keeping all cells.")
    seu$miQC_prob <- NA_real_
    seu$keep_miQC <- TRUE
    return(seu)
  }

  # Fit miQC (linear → spline → 1D)
  set.seed(1234)
  mmod <- try(miQC::mixtureModel(sce, model_type = "linear"), silent = TRUE)
  if (inherits(mmod, "try-error")) {
    set.seed(5678)
    mmod <- try(miQC::mixtureModel(sce, model_type = "spline"), silent = TRUE)
  }
  if (inherits(mmod, "try-error")) {
    set.seed(91011)
    mmod <- try(miQC::mixtureModel(sce, model_type = "one_dimensional"), silent = TRUE)
  }
  if (inherits(mmod, "try-error")) {
    warning("miQC::mixtureModel failed after all model types; keeping all cells.")
    seu$miQC_prob <- NA_real_
    seu$keep_miQC <- TRUE
    return(seu)
  }

  # Posterior & choose “compromised” component (higher mito at median detected)
  post_mat <- flexmix::posterior(mmod)
  coef_mat <- flexmix::parameters(mmod)
  det_med <- stats::median(SummarizedExperiment::colData(sce)$detected, na.rm = TRUE)
  pick_comp <- function(coefs, det_med) {
    if (is.null(dim(coefs))) {
      return(2L)
    } # one_dimensional
    if (all(c("(Intercept)", "detected") %in% rownames(coefs))) {
      pred <- coefs["(Intercept)", ] + det_med * coefs["detected", ]
      return(which.max(pred)) # higher predicted mito% = compromised
    } else {
      return(2L)
    }
  }
  comp_k <- pick_comp(coef_mat, det_med)

  post_comp <- post_mat[, comp_k, drop = TRUE]
  names(post_comp) <- if (!is.null(rownames(post_mat))) rownames(post_mat) else colnames(sce)

  # Align back to Seurat barcodes
  seu$miQC_prob <- post_comp[match(colnames(seu), names(post_comp))]

  # Keep criterion: keep cells with posterior < 0.75 (miQC default)
  seu$keep_miQC <- seu$miQC_prob < 0.75

  if (isTRUE(debug)) {
    dir.create("qc_outputs", showWarnings = FALSE, recursive = TRUE)
    df <- data.frame(
      barcode = colnames(sce),
      detected = SummarizedExperiment::colData(sce)$detected,
      subsets_mito_percent = SummarizedExperiment::colData(sce)$subsets_mito_percent,
      miQC_prob = seu$miQC_prob[match(colnames(sce), colnames(seu))],
      keep_miQC = seu$keep_miQC[match(colnames(sce), colnames(seu))],
      stringsAsFactors = FALSE
    )
    utils::write.csv(df, file = "qc_outputs/miqc_input_metrics.csv", row.names = FALSE)
  }
  return(seu)
}


# ===========================
# 7) Baseline: DropletQC damaged-only (requires nuclear_fraction + umi)
# ===========================

baseline_dropletqc <- function(seu, seurat_sp) {
  ## 1) Grab counts for spliced/unspliced and compute per-cell totals
  spliced <- GetAssayData(seurat_sp, assay = "spliced", layer = "counts")
  unspliced <- GetAssayData(seurat_sp, assay = "unspliced", layer = "counts")

  # Sanity: same cells in same order
  stopifnot(identical(colnames(spliced), colnames(unspliced)))

  exon_sum <- Matrix::colSums(spliced) # exonic counts
  intron_sum <- Matrix::colSums(unspliced) # intronic counts

  # DropletQC nuclear fraction = intronic / (intronic + exonic)
  nf <- intron_sum / (intron_sum + exon_sum + 1e-12)
  umi <- intron_sum + exon_sum # UMI for nf/umi plane (do NOT include ambiguous here)

  ## 2) Add to Seurat metadata
  seurat_sp$nuclear_fraction <- as.numeric(nf[colnames(seurat_sp)])
  seurat_sp$umi_nf_plane <- as.numeric(umi[colnames(seurat_sp)])

  # 1) Build the data frame for identify_empty_drops
  nf_umi <- data.frame(
    nf = as.numeric(seurat_sp$nuclear_fraction),
    umi = as.numeric(seurat_sp$umi_nf_plane),
    row.names = colnames(seurat_sp)
  )

  # Optional: drop any rows with NA (NF/UMI)
  nf_umi <- nf_umi[complete.cases(nf_umi[, 1:2]), ]

  # 2) Run identify_empty_drops
  ed <- identify_empty_drops(nf_umi = nf_umi, include_plot = FALSE)

  # 3) Add cell_type for identify_damaged_cells using the standard object's clusters
  # Use .normalize_cells internally if names dont perfectly match, but here they should due to subset alignment
  ed$cell_type <- seu$seurat_clusters[match(.normalize_cells(rownames(ed)), .normalize_cells(colnames(seu)))]

  # Drop any cell types that mapped to NA before pushing into identify_damaged_cells, else EM will crash
  ed <- ed[!is.na(ed$cell_type), ]

  # 4) Run DropletQC’s damaged-cell classifier
  dc <- identify_damaged_cells(
    ed, # min NF separation between intact vs damaged (tune if needed)
    verbose = TRUE,
    output_plots = T
  )
  dc$plots
  table(dc[[1]]$cell_status)

  # 5) Save calls back to Seurat metadata
  calls <- dc[[1]] # first list element = annotated table
  seurat_sp$DropletQC_status <- calls$cell_status[match(colnames(seurat_sp), rownames(calls))]
  seurat_sp$DropletQC_status[is.na(seurat_sp$DropletQC_status)] <- "unknown_NA"

  # labels can be "cell", "empty_droplet", or "damaged_cell"
  kept_dq <- seurat_sp$DropletQC_status == "cell"
  kept_dq[is.na(kept_dq)] <- FALSE
  seurat_sp$keep_dropletqc <- kept_dq

  want <- c("nuclear_fraction", "umi_nf_plane", "DropletQC_status", "keep_dropletqc")
  seu_merged_meta <- merge_meta_only(seu, seurat_sp, cols = want, suffix = "")
  return(seu_merged_meta)
}


## qc plot -----------------------------------------------------------------


seurat_std_small <- baseline_seurat(seurat_std_small)
table(seurat_std_small$keep_seurat)

p <- DimPlot(seurat_std_small, group.by = "keep_seurat", cols = c("red", "gray"), pt.size = 0.5, alpha = 0.5) +
  ggtitle("Seurat") +
  NoLegend()
ggsave("keep_seurat.png", p, width = 4, height = 4)

seurat_std_small <- baseline_scater(seurat_std_small)
table(seurat_std_small$keep_scater)

p <- DimPlot(seurat_std_small, group.by = "keep_scater", cols = c("red", "gray"), pt.size = 0.5, alpha = 0.5) +
  ggtitle("scater") +
  NoLegend()
ggsave("keep_scater.png", p, width = 4, height = 4)


seurat_std_small <- baseline_miQC(seurat_std_small)
if (table(seurat_std_small$keep_miQC)[1] > table(seurat_std_small$keep_miQC)[2]) {
  seurat_std_small$keep_miQC <- !seurat_std_small$keep_miQC
}
table(seurat_std_small$keep_miQC)

p <- DimPlot(seurat_std_small, group.by = "keep_miQC", cols = c("red", "gray"), pt.size = 0.5, alpha = 0.5) +
  ggtitle("miQC") +
  NoLegend()
p
ggsave("keep_miqc.png", p, width = 4, height = 4)
table(seurat_std_small$keep_miQC)


seurat_std_small <- baseline_dropletqc(seurat_std_small, seurat_sp_small)
table(seurat_std_small$keep_dropletqc)
p <- DimPlot(seurat_std_small, group.by = "keep_dropletqc", cols = c("red", "gray"), pt.size = 0.5, alpha = 0.5) +
  ggtitle("DropletQC") +
  NoLegend()
p
ggsave("keep_dq.png", p, width = 4, height = 4)


# run scqc ----------------------------------------------------------------


library(scQCenrich)
res <- run_qc_pipeline(
  obj = seurat_std_small, # normal Seurat for baseline QC
  species = organism,
  assay = "RNA",
  method = "gmm",
  # external splicing sources:
  spliced_obj = seurat_sp_small,
  spliced_assay = "spliced",
  spliced_layer = "counts",
  unspliced_obj = seurat_sp_small,
  unspliced_assay = "unspliced",
  unspliced_layer = "counts",
  report_html = T,
  report_file = "qc_outputs/qc_report.html",
  debug = TRUE,
  annot_method = "none",
  tissue = organ,
  marker_method = "findmarkers",
  doublets = "none",
  enrichment_plots = T
  # qc_strength = "strict"
)
if (!is.null(res$report)) browseURL(res$report)

seu <- res$obj_all
seu$keep_enrich <- seu$qc_status != "remove"
seurat_std_small$keep_enrich <- seu$keep_enrich[match(colnames(seurat_std_small), colnames(seu))]


p <- DimPlot(seu, group.by = "keep_enrich", cols = c("red", "gray"), pt.size = 0.5, alpha = 0.5) +
  ggtitle("scQCenrich") +
  NoLegend()
p
ggsave("keep_enrich.png", p, width = 4, height = 4)




# ---- scQCenrich WITHOUT splice data -------------------------


cat("\n--- Running scQCenrich without splice data (MALAT1 proxy) ---\n")
met_nosplice <- calcQCmetrics(
  obj         = seurat_std_small,
  species     = organism,
  assay       = "RNA",
  add_to_meta = FALSE
  # spliced_obj intentionally omitted -> MALAT1 proxy activates
)
cat(
  "intronic_frac all-NA (confirming no splice used):",
  all(is.na(met_nosplice$intronic_frac)), "\n"
)
cat("MALAT1_frac range:", round(range(met_nosplice$MALAT1_frac, na.rm = TRUE), 4), "\n")

res_nosplice <- flagLowQuality(
  metrics     = met_nosplice,
  method      = "gmm",
  qc_strength = "auto",
  rescue_mode = "lenient"
)

# Align barcodes back to seurat_std_small (metrics rownames == colnames(seu))
seurat_std_small$keep_enrich_nosplice <- res_nosplice$qc_status[
  match(colnames(seurat_std_small), rownames(res_nosplice))
] != "remove"

cat("kept (no splice):", sum(seurat_std_small$keep_enrich_nosplice, na.rm = TRUE), "\n")
cat("removed (no splice):", sum(!seurat_std_small$keep_enrich_nosplice, na.rm = TRUE), "\n")

p_nosplice <- DimPlot(
  seurat_std_small,
  group.by = "keep_enrich_nosplice",
  cols = c("red", "gray"), pt.size = 0.5, alpha = 0.5
) +
  ggtitle("scQCenrich (no splice)") +
  NoLegend()
p_nosplice
ggsave("keep_enrich_nosplice.png", p_nosplice, width = 4, height = 4)
