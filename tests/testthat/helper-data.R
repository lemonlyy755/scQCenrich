load_real_objs <- function() {
  std_path <- testthat::test_path("testdata", "seurat_std_small.rds")
  sp_path  <- testthat::test_path("testdata", "seurat_sp_small.rds")

  message("[scQCenrich tests] Looking for data in: ", testthat::test_path("testdata"))
  message("  - std: ", std_path, " (exists=", file.exists(std_path), ")")
  message("  -  sp: ", sp_path,  " (exists=", file.exists(sp_path),  ")")

  # fallback: also accept old location tests/testdata and inst/testdata
  if (!file.exists(std_path)) {
    alt1 <- file.path("tests", "testdata", "seurat_std_small.rds")
    alt2 <- system.file("testdata", "seurat_std_small.rds", package = "scQCenrich")
    for (p in c(alt1, alt2)) if (nzchar(p) && file.exists(p)) { std_path <- p; break }
  }
  if (!file.exists(sp_path)) {
    alt1 <- file.path("tests", "testdata", "seurat_sp_small.rds")
    alt2 <- system.file("testdata", "seurat_sp_small.rds", package = "scQCenrich")
    for (p in c(alt1, alt2)) if (nzchar(p) && file.exists(p)) { sp_path <- p; break }
  }

  testthat::skip_if_not(file.exists(std_path), "seurat_std_small.rds missing")
  testthat::skip_if_not(file.exists(sp_path),  "seurat_sp_small.rds missing")

  std <- readRDS(std_path)
  sp  <- readRDS(sp_path)

  stopifnot("RNA" %in% Seurat::Assays(std))
  stopifnot(all(c("spliced","unspliced") %in% Seurat::Assays(sp)))

  list(std = std, sp = sp)
}
