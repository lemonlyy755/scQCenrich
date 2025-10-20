

# man ---------------------------------------------------------------------

usethis::use_logo("man/figures/logo.png")

usethis::use_git_ignore("dev.R")


# combine script ----------------------------------------------------------


# Combine all R scripts in R/ into one text file
combine_package_files <- function(output_file = "package_functions.txt") {
  # List all .R files in R/ folder
  files <- list.files("R", full.names = TRUE)
  files <- append(files,
                  list.files("inst/templates", full.names = TRUE)
  )

  # Read them all and add headers for clarity
  combined <- unlist(lapply(files, function(f) {
    c(
      paste0("# ---- File: ", basename(f), " ----"),
      readLines(f),
      ""
    )
  }))

  # Write to a single text file
  writeLines(combined, output_file)
  message("Combined into: ", output_file)
}

# Run
combine_package_files()


#  combine all ------------------------------------------------------------

# ---- File: export_pkg_txt.R ----
# Combine all package files into a single text file, skipping specified folders or file types

combine_pkg_files <- function(pkg_dir = ".",
                              out_file = "package_combined.txt",
                              skip_files = NULL,
                              skip_dirs = c(".git", ".github", "data", "man", "tests", "docs"),
                              skip_ext  = c(".rda", ".rds", ".png", ".jpg", ".jpeg", ".pdf", ".html", ".csv", ".tsv")) {

  # --- find all files recursively ---
  files <- list.files(pkg_dir, recursive = TRUE, full.names = TRUE)

  # --- skip folders ---
  if (length(skip_dirs) > 0) {
    skip_pattern <- paste0("(^|/)", skip_dirs, "(/|$)", collapse = "|")
    files <- files[!grepl(skip_pattern, files)]
  }

  # --- skip by extension ---
  if (length(skip_ext) > 0) {
    ext_pattern <- paste0("(", paste0("\\", skip_ext, collapse = "|"), ")$")
    files <- files[!grepl(ext_pattern, files, ignore.case = TRUE)]
  }

  # --- skip individual files ---
  if (length(skip_files) > 0) {
    files <- files[!basename(files) %in% skip_files]
  }

  cat("Found", length(files), "files to include\n")

  # --- write combined file ---
  out_con <- file(out_file, "w", encoding = "UTF-8")
  for (f in files) {
    cat("\n# ---- File:", f, "----\n", file = out_con, append = TRUE)
    lines <- readLines(f, warn = FALSE, encoding = "UTF-8")
    writeLines(lines, out_con)
    cat("\n", file = out_con, append = TRUE)
  }
  close(out_con)
  message("Written combined file to: ", out_file)
}

# Example run:
combine_pkg_files(
  skip_files = c(
    "dev.R", "toy_seu.rda", "README.md",
    "testloom.rds", "seurat_std.rds", "PanglaoDB_markers_27_Mar_2020.tsv",
    "package_combined.txt"
  ),
  skip_dirs = c(".git", "data", "man", "inst/extdata", "docs", ".github","qc_outputs","benchmarking"),
  skip_ext  = c(".rda", ".rds", ".png", ".jpg", ".jpeg", ".pdf", ".html", ".csv", ".tsv",".log",".txt",".Rhistory",".Rproj")
)

# run script --------------------------------------------------------------
options(timeout = 600)
options(scQCenrich.debug = TRUE, scQCenrich.debug_file = "qc_outputs/sctype_debug.log")

library(devtools)
options(repos = BiocManager::repositories())

document()
install()





library(scQCenrich)
res <- run_qc_pipeline(
  obj          = seurat_std_small,                 # normal Seurat for baseline QC
  species      = "mouse",
  assay        = "RNA",
  method       = "gmm",
  spliced_obj        = seurat_sp_small,
  spliced_assay      = "spliced",
  spliced_layer      = "counts",
  unspliced_obj      = seurat_sp_small,
  unspliced_assay    = "unspliced",
  unspliced_layer    = "counts",
  report_html  = TRUE,
  report_file  = "qc_outputs/qc_report.html",
  debug        = TRUE,
  annot_method = "marker_score",
  tissue = c("Skin","Immune system","Connective tissue"),
  marker_method = "findmarkers",
  doublets  = "none",
  enrichment_plots = F
  #qc_strength = "strict"


)

 browseURL(res$report)
 head(res$obj_all@meta.data)
 table(res$status_df$qc_status)
 # expect non-zero counts in 'borderline' and/or 'remove'
 head(subset(res$status_df, qc_status != "keep"))


 traceback()

 devtools::check()


# finish up ---------------------------------------------------------------

 # roxygen docs & NAMESPACE
 devtools::document()
 devtools::load_all()

 devtools::test()
 # build vignettes into inst/doc/ and run checks
 devtools::build(vignettes = TRUE)
 devtools::install(build_vignettes = TRUE)
 devtools::check()


 devtools::build_vignettes()


 files <- list.files("R", pattern = "\\.[Rr]$", full.names = TRUE)
 lapply(files, tools::showNonASCIIfile)

# test subset helper ------------------------------------------------------

 library(Seurat)
 seurat_std <- readRDS("seurat_std.rds")
 seurat_sp <- readRDS("testloom.rds")
 organism <- "mouse"
 organ <- c("Skin","Immune system","Connective tissue")


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
 sp_ids  <- .normalize_cells(colnames(seurat_sp))

 # Find intersection
 common <- intersect(std_ids, sp_ids)
 message("Found ", length(common), " overlapping cells.")

 # Sample 500 cells (or fewer if less available)
 set.seed(123)  # for reproducibility
 keep_ids <- sample(common, min(500, length(common)))
 keep_ids <- sample(common, min(10000, length(common)))

 # Map back to original names in each object
 std_keep <- colnames(seurat_std)[.normalize_cells(colnames(seurat_std)) %in% keep_ids]
 sp_keep  <- colnames(seurat_sp)[ .normalize_cells(colnames(seurat_sp))  %in% keep_ids]

 # Subset
 seurat_std_small <- subset(seurat_std, cells = std_keep)
 seurat_sp_small  <- subset(seurat_sp,  cells = sp_keep)

 # Check dimensions
 dim(seurat_std_small)  # genes x 500 (or fewer)
 dim(seurat_sp_small)

# auto generation ---------------------------------------------------------




 usethis::use_vignette("scQCenrich-intro")

 devtools::build_readme()  # knits README.Rmd -> README.md

# create toy data ---------------------------------------------------------




 suppressPackageStartupMessages({
   library(org.Mm.eg.db) # or org.Hs.eg.db
   library(AnnotationDbi)
   library(Matrix); library(Seurat); library(usethis)
 })

 set.seed(1)
 ng <- 1000L; nc <- 500L

 # sample known SYMBOLs (mouse shown here)
 sym <- AnnotationDbi::keys(org.Mm.eg.db, keytype = "SYMBOL")
 genes <- sample(sym, ng, replace = FALSE)

 mat_dense <- matrix(rpois(ng * nc, lambda = 0.2), nrow = ng, ncol = nc,
                     dimnames = list(genes, paste0("cell", seq_len(nc))))
 counts <- Matrix::Matrix(mat_dense, sparse = TRUE)

 toy_seu <- Seurat::CreateSeuratObject(counts)
 toy_seu <- Seurat::NormalizeData(toy_seu, verbose = FALSE)
 usethis::use_data(toy_seu, overwrite = TRUE, compress = "xz")




# marker annotaion debug --------------------------------------------------



 options(scQCenrich.debug = TRUE)

 # choose one: generic signatures OR panglao_signatures_debug(...)
 sig <- NULL
  sig <- panglao_signatures_debug(tsv = "inst/extdata/PanglaoDB_markers_27_Mar_2020.tsv", species = "mouse", ui_max = 0.20,
                                  canonical_only = TRUE, min_genes = 2,
                                  obj = seurat_std_small, assay = "RNA")

 res <- .annot_marker_score_debug(
   obj = seurat_std_small,
   species = "mouse",
   assay   = "RNA",
   level   = "cluster",            # or "cell"
   signatures = sig,               # NULL = use generic_signatures()
   unknown_min_genes  = 2,
   unknown_min_margin = 0.15,
   unknown_min_top    = 0.10,
   nbin = 24, ctrl = 100,          # start from Seurat defaults; tweak later
   autotune_thresholds = TRUE,
   canonical_only = F,
   min_genes = 2
 )

# git ---------------------------------------------------------------------
 # in a clean working tree
 git checkout --orphan fresh-start
 git add -A
 git commit -m "Fresh start"
 git branch -M main
 git remote remove origin 2>/dev/null || true
 git remote add origin https://github.com/lemonlyy755/scQCenrich.git
 git push -u --force origin main


 # 2) Stage all tracked + new files
 git add -A

 # 3) Commit with a clear message
 git commit -m "pass check 3"



# comaprison --------------------------------------------------------------

 library(Seurat)
 pbmc_filt <- Read10X_h5("filtered_feature_bc_matrix.h5") |> CreateSeuratObject()




