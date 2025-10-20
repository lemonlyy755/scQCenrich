# data-raw/toy_seu.R
# Create a small synthetic Seurat object (~500 cells, 1000 genes) and save as data/toy_seu.rda

suppressPackageStartupMessages({
  library(Matrix)   # for sparse matrices
  library(Seurat)   # for CreateSeuratObject/NormalizeData
  library(usethis)  # for use_data()
})

set.seed(1)

ng <- 1000L   # genes
nc <- 500L    # cells

# Small, sparse counts with low lambda to mimic single-cell sparsity
mat_dense <- matrix(rpois(ng * nc, lambda = 0.2), nrow = ng, ncol = nc)
rownames(mat_dense) <- paste0("Gene", seq_len(ng))
colnames(mat_dense) <- paste0("cell", seq_len(nc))
counts <- Matrix::Matrix(mat_dense, sparse = TRUE)

# Build Seurat obj + quick normalization so downstream examples work
toy_seu <- Seurat::CreateSeuratObject(counts)
toy_seu <- Seurat::NormalizeData(toy_seu, verbose = FALSE)

# Save into the package's data/ folder
usethis::use_data(toy_seu, overwrite = TRUE, compress = "xz")
