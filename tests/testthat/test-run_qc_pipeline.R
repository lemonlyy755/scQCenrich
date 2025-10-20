test_that("pipeline runs on real data (minimal)", {
  objs <- load_real_objs()
  std <- objs$std
  sp  <- objs$sp

  expect_no_error({
    res <- run_qc_pipeline(
      obj       = std,
      species   = "mouse",
      assay     = "RNA",
      method    = "gmm",
      # external splicing inputs:
      spliced_obj     = sp,
      spliced_assay   = "spliced",
      spliced_layer   = "counts",
      unspliced_obj   = sp,
      unspliced_assay = "unspliced",
      unspliced_layer = "counts",
      report_html = FALSE,
      debug       = TRUE,
      annot_method = "marker_score",
      tissue       = "Skin",
      doublets     = "none"
    )
  })
})
