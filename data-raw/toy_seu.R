# data-raw/toy_seu.R
# Create a small PBMC-derived Seurat object (~500 cells, ~13k genes)
# from benchmark data and save as data/toy_seu.rda

suppressPackageStartupMessages({
  library(Matrix)
  library(Seurat)
  library(usethis)
})

set.seed(1)

src_path <- file.path("benchmarking", "qc_pbmc", "pbmc_annotated.rds")
stopifnot(file.exists(src_path))

src <- readRDS(src_path)
stopifnot(
  inherits(src, "Seurat"),
  "cell_type" %in% colnames(src@meta.data)
)

# Sample a balanced subset from the 10 largest benchmark PBMC cell types.
celltype_tab <- sort(table(src$cell_type), decreasing = TRUE)
keep_types <- names(celltype_tab)[seq_len(min(10L, length(celltype_tab)))]
cells_use <- unlist(lapply(keep_types, function(tp) {
  avail <- colnames(src)[src$cell_type == tp]
  sample(avail, min(50L, length(avail)))
}), use.names = FALSE)

toy_seu <- subset(src, cells = cells_use)
counts <- Seurat::GetAssayData(toy_seu, assay = "RNA", layer = "counts")
keep_genes <- rownames(counts)[Matrix::rowSums(counts > 0) >= 5]

# Keep only the RNA assay plus counts/data layers to keep the bundled dataset compact.
toy_seu <- Seurat::DietSeurat(
  toy_seu,
  assays = "RNA",
  layers = c("counts", "data"),
  features = keep_genes,
  dimreducs = NULL,
  graphs = NULL,
  misc = FALSE
)

toy_seu@meta.data <- toy_seu@meta.data[
  ,
  intersect(c("orig.ident", "cell_type"), colnames(toy_seu@meta.data)),
  drop = FALSE
]

# Save into the package's data/ folder
usethis::use_data(toy_seu, overwrite = TRUE, compress = "xz")
