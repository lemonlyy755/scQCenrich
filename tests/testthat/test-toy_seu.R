test_that("bundled toy_seu supports marker-based annotation", {
  data("toy_seu", package = "scQCenrich")

  expect_s4_class(toy_seu, "Seurat")
  expect_equal(ncol(toy_seu), 500L)
  expect_true("cell_type" %in% colnames(toy_seu@meta.data))

  expect_no_error({
    annotated <- suppressWarnings(auto_annotate(
      toy_seu,
      species = "human",
      assay = "RNA",
      prefer = "marker_score",
      level = "cluster",
      marker_score_level = "cluster"
    ))
  })

  expect_true("auto_celltype" %in% colnames(annotated@meta.data))
  expect_gt(sum(!is.na(annotated$auto_celltype)), 0)
})
