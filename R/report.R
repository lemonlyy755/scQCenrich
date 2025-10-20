# ---- R/report.R ----
# Render a pretty HTML QC report from an external Rmd template
# If `template` is NULL, we copy the bundled template to `save_dir/qc_report.Rmd`.

#' Render a pretty HTML QC report
#' @param obj Seurat object (with qc_status in meta)
#' @param metrics data.frame from calcQCmetrics()
#' @param status_df data.frame from flagLowQuality()
#' @param save_dir output dir for assets/plots (and where the Rmd will be copied)
#' @param report_file output HTML path (default: file.path(save_dir, "qc_report.html"))
#' @param theme_bootswatch Bootswatch theme name
#' @param self_contained embed all assets into single HTML (TRUE)
#' @param sample_col optional meta column for per-sample summaries
#' @param template path to a .Rmd template; if NULL, copy bundled template to save_dir and use it
#' @param species Character scalar, one of "mouse" or "human"
#' @export
render_qc_report <- function(
    obj, metrics, status_df,
    save_dir    = "qc_outputs",
    report_file = file.path(save_dir, "qc_report.html"),
    theme_bootswatch = "minty",
    self_contained = TRUE,
    sample_col = NULL,
    template = NULL,
    species = c("mouse","human")
){
  species <- match.arg(species)
  if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE, showWarnings = TRUE)
  out_dir <- normalizePath(dirname(report_file), winslash = "/", mustWork = FALSE)
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  have_rmd <- requireNamespace("rmarkdown", quietly = TRUE)

  if (is.null(template)) {
    template <- .ensure_qc_report_template(save_dir)
  } else {
    if (!file.exists(template)) stop("Template not found: ", template)
  }

  if (!have_rmd) {
    warning("rmarkdown not available; writing a minimal HTML via htmltools.")
    html <- htmltools::tagList(
      htmltools::h2("scQCenrich QC Report"),
      htmltools::p(sprintf("Cells kept: %d; removed: %d",
                           sum(status_df$qc_status == "keep", na.rm = TRUE),
                           sum(status_df$qc_status == "remove", na.rm = TRUE)))
    )
    htmltools::save_html(html, file = report_file, background = "white")
    return(invisible(report_file))
  }

  env <- new.env(parent = globalenv())
  assign("obj",        obj,        envir = env)
  assign("metrics",    metrics,    envir = env)
  assign("status_df",  status_df,  envir = env)
  assign("save_dir",   normalizePath(save_dir, winslash = "/", mustWork = FALSE), envir = env)
  assign("SCOL",       sample_col, envir = env)
  assign("theme_bootswatch", theme_bootswatch, envir = env)
  assign("species",    species,    envir = env)

  ofile <- basename(report_file)
  rmarkdown::render(
    input         = template,
    output_file   = ofile,
    output_dir    = out_dir,
    output_format = rmarkdown::html_document(self_contained = self_contained),
    envir         = env,
    quiet         = FALSE
  )
  invisible(file.path(out_dir, ofile))
}

.ensure_qc_report_template <- function(save_dir) {
  # try installed path first (inst/templates/qc_report.Rmd)
  pkg <- try(utils::packageName(), silent = TRUE)
  bundled <- if (!inherits(pkg, "try-error") && nzchar(pkg)) {
    system.file("templates", "qc_report.Rmd", package = pkg)
  } else ""
  if (!nzchar(bundled) || !file.exists(bundled)) {
    candidate <- file.path("inst", "templates", "qc_report.Rmd")
    if (file.exists(candidate)) bundled <- candidate
  }
  if (!file.exists(bundled)) stop("Bundled qc_report.Rmd not found. Did you add it to inst/templates/?")

  target <- file.path(save_dir, "qc_report.Rmd")
  if (!file.exists(target)) {
    dir.create(save_dir, showWarnings = FALSE, recursive = TRUE)
    file.copy(bundled, target, overwrite = FALSE)
    message("Copied qc_report.Rmd to: ", normalizePath(target, winslash = "/"))
  }
  normalizePath(target, winslash = "/", mustWork = FALSE)
}
