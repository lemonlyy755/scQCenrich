# ---- R/utils.R ----
# Utilities + robust Seurat v5 layer handling and debugging

`%||%` <- function(x, y) if (is.null(x)) y else x  # small local helper

.dbg <- function(...) {
  if (isTRUE(getOption("scQCenrich.debug", FALSE))) {
    message("[scQCenrich-debug] ", paste(..., collapse = ""))
  }
}

.mad_thresholds <- function(x, k = 3) {
  m  <- stats::median(x, na.rm = TRUE)
  md <- stats::mad(x, constant = 1, na.rm = TRUE)
  list(
    high = m + k * md,
    low = max(0, m - k * md),
    med = m,
    mad = md
  )
}

.to_dgC <- function(x, name = "unknown") {
  if (inherits(x, "dgCMatrix")) {
    .dbg(name, ": already dgCMatrix (", paste(dim(x), collapse = "x"), ")")
    return(x)
  }
  out <- suppressWarnings(try(methods::as(x, "dgCMatrix"), silent = TRUE)
  )
  if (!inherits(out, "try-error") && inherits(out, "dgCMatrix"))
  {
    .dbg(name,
         ": coerced via methods::as -> dgCMatrix (",
         paste(dim(out), collapse = "x"),
         ")")
    return(out)
  }
  m <- suppressWarnings(try(as.matrix(x), silent = TRUE)
  )
  if (!inherits(m, "try-error") && is.matrix(m))
  {
    out <- Matrix::Matrix(m, sparse = TRUE)
    out <- methods::as(out, "dgCMatrix")
    .dbg(
      name,
      ": coerced via as.matrix -> sparse dgCMatrix (",
      paste(dim(out), collapse = "x"),
      ")"
    )
    return(out)
  }
  if (is.vector(x) && is.atomic(x)) {
    m <- matrix(x, ncol = 1)
    out <- Matrix::Matrix(m, sparse = TRUE)
    out <- methods::as(out, "dgCMatrix")
    .dbg(name,
         ": vector->matrix -> dgCMatrix (",
         paste(dim(out), collapse = "x"),
         ")")
    return(out)
  }
  stop(
    "Cannot coerce object of class [",
    paste(class(x), collapse = ","),
    "] to dgCMatrix (",
    name,
    ")."
  )
}

.parse_assay_layer <- function(x, default_layer = "counts") {
  if (is.null(x))
    return(list(assay = NULL, layer = NULL))
  if (length(x) != 1L)
    stop(".parse_assay_layer(): expecting a single string.")
  if (grepl(":", x, fixed = TRUE)) {
    parts <- strsplit(x, ":", fixed = TRUE)[[1]]
    if (length(parts) != 2L)
      stop("Use 'assay:layer', e.g., 'RNA:unspliced'.")
    list(assay = parts[1], layer = parts[2])
  } else {
    list(assay = x, layer = default_layer)
  }
}

# v5/v4-safe getter that ALWAYS returns a dgCMatrix
.get_assay_data <- function(obj, assay, layer_or_slot) {
  .dbg("GetAssayData(", assay, " : ", layer_or_slot, ")")
  # 1) Preferred v5: explicit layer on assay
  out <- suppressWarnings(try(SeuratObject::GetAssayData(obj, assay = assay, layer = layer_or_slot),
                              silent = TRUE)
  )
  if (!inherits(out, "try-error") && !is.null(out))
  {
    return(.to_dgC(out, paste0(
      assay, ":", layer_or_slot, " [v5-layer]"
    )))
  }
  # 2) Try using assay object directly
  ao <- suppressWarnings(try(obj[[assay]], silent = TRUE)
  )
  if (!inherits(ao, "try-error") && !is.null(ao))
  {
    out <- suppressWarnings(try(SeuratObject::GetAssayData(ao, layer = layer_or_slot),
                                silent = TRUE)
    )
    if (!inherits(out, "try-error") && !is.null(out))
    {
      return(
        .to_dgC(out, paste0(
          assay, ":", layer_or_slot, " [assay-layer]"
        )
        ))
    }
  }
  # 3) v4 fallback: slot (deprecated in v5)
  out <- suppressWarnings(
    SeuratObject::GetAssayData(obj, assay = assay, slot = layer_or_slot)
  )
  .to_dgC(
    out, paste0(assay, ":", layer_or_slot, " [slot-fallback]")
  )
}

#' @keywords internal
.ensure_layers_and_data <- function(obj,
                                    assay = "RNA",
                                    debug_dir = NULL) {
  Seurat::DefaultAssay(obj) <- assay
  if (utils::packageVersion("SeuratObject") >= "5.0.0") {
    obj <- tryCatch(
      SeuratObject::JoinLayers(obj, assay = assay),
      error = function(e) {
        if (!is.null(debug_dir))
          writeLines(e$message,
                     file.path(debug_dir, "joinlayers_error.txt"))
        obj
      }
    )
  }
  dat <- suppressWarnings(try(SeuratObject::GetAssayData(obj, assay = assay, layer = "data"),
                              silent = TRUE)
  )
  if (inherits(dat, "try-error") ||
      is.null(dat) || nrow(dat) == 0L || ncol(dat) == 0L)
  {
    obj <- Seurat::NormalizeData(obj, verbose = FALSE)
  }
  if (!is.null(debug_dir)) {
    dir.create(debug_dir, FALSE, TRUE)
    info <- list(
      R.version       = R.version.string,
      Seurat          = as.character(utils::packageVersion("Seurat")),
      SeuratObject    = as.character(utils::packageVersion("SeuratObject")),
      assay           = assay,
      n_features_cnts = nrow(Seurat::GetAssayData(
        obj, assay = assay, layer = "counts"
      )),
      n_cells_cnts    = ncol(Seurat::GetAssayData(
        obj, assay = assay, layer = "counts"
      )),
      n_features_data = nrow(Seurat::GetAssayData(
        obj, assay = assay, layer = "data"
      )),
      n_cells_data    = ncol(Seurat::GetAssayData(
        obj, assay = assay, layer = "data"
      ))
    )
    write.csv(
      as.data.frame(info, check.names = FALSE),
      file.path(debug_dir, "00_assay_layer_summary.csv"),
      row.names = FALSE
    )
  }
  obj
}

#' @keywords internal
.log_cpm <- function(m,
                     prior.count = 1,
                     scale = 1e6) {
  m <- .to_dgC(m, "log_cpm input")
  if (nrow(m) == 0L || ncol(m) == 0L) {
    .dbg("log_cpm: empty matrix -> returning empty matrix")
    return(m)
  }
  libs <- Matrix::colSums(m)
  .dbg(sprintf("log_cpm: lib length=%d", length(libs)))
  libs[!is.finite(libs)] <- 0

  # preallocate zeros; keep dimnames
  out <- Matrix::Matrix(
    0,
    nrow = nrow(m),
    ncol = ncol(m),
    sparse = TRUE,
    dimnames = dimnames(m)
  )

  keep <- libs > 0
  if (!any(keep)) {
    .dbg("log_cpm: all libraries are zero -> returning all zeros")
    return(out)
  }

  # log1p((counts + prior)/lib * scale), only for columns with nonzero libs
  S <- Matrix::Diagonal(x = scale / libs[keep])
  out[, keep] <- log1p((m[, keep, drop = FALSE] + prior.count) %*% S)
  out
}

.safe_ggsave <- function(path,
                         plot,
                         width = 6.5,
                         height = 5,
                         dpi = 200) {
  try({
    if (!is.null(plot))
      ggplot2::ggsave(path,
                      plot,
                      width = width,
                      height = height,
                      dpi = dpi)
  }, silent = TRUE)
  invisible(path)
}

.guess_id_type <- function(genes) {
  g <- genes[!is.na(genes)]
  frac_hs  <- mean(grepl("^ENSG\\d+(?:\\.\\d+)?$", g))
  frac_mm  <- mean(grepl("^ENSMUSG\\d+(?:\\.\\d+)?$", g))
  if (max(frac_hs, frac_mm, na.rm = TRUE) > 0.6)
    "ENSEMBL"
  else
    "SYMBOL"
}

.strip_ver <- function(x)
  sub("\\.\\d+$", "", x)

#' @export
# internal helper; DO NOT export (or keep your current export status)
.mt_gene_indices <- function(genes, species = c("mouse", "human")) {
  species <- match.arg(species)

  # Fast path for SYMBOL rownames: match MT- or mt.  (e.g., "MT-CO1")
  i <- grep("^mt[-\\.]", genes, ignore.case = TRUE)
  if (length(i)) return(i)

  # Only attempt DB mapping if rownames *look* like Ensembl
  looks_ens <- any(grepl("^(ENSG\\d+|ENSMUSG\\d+)", genes))
  if (!looks_ens) return(integer(0))

  # Graceful fallbacks if DBs are missing (avoid errors in light tests)
  if (!requireNamespace("AnnotationDbi", quietly = TRUE)) return(integer(0))
  pkg <- if (species == "mouse") "org.Mm.eg.db" else "org.Hs.eg.db"
  if (!requireNamespace(pkg, quietly = TRUE)) return(integer(0))

  # Get the OrgDb object from the namespace safely
  db <- getExportedValue(pkg, pkg)

  # Map Ensembl -> SYMBOL and then detect MT- symbols
  suppressWarnings({
    map <- AnnotationDbi::select(
      db,
      keys    = genes,
      keytype = "ENSEMBL",
      columns = "SYMBOL"
    )
  })
  if (!NROW(map)) return(integer(0))

  sym <- setNames(map$SYMBOL, map$ENSEMBL)
  sym <- sym[match(genes, names(sym))]
  which(grepl("^MT[-\\.]", sym, ignore.case = TRUE))
}


# canonicalize common barcode variants
.norm_cell_id <- function(x) {
  x <- as.character(x)
  x <- sub("^.*:", "", x)                  # drop "sample:" prefix
  x <- sub("[\\-_.:][0-9A-Za-z]+$", "", x) # drop -1/_1/:1/etc.
  x <- sub("x$", "", x, ignore.case = TRUE)
  x <- toupper(x)
  gsub("[^ACGT]", "", x)
}

# build mapping from base cells -> external by normalized barcode
.align_external_by_cell <- function(base_cells,
                                    ext_cells,
                                    debug = getOption("scQCenrich.debug", FALSE)) {
  bkey <- .norm_cell_id(base_cells)
  ekey <- .norm_cell_id(ext_cells)

  # map base keys to FIRST unique external index; ambiguous external keys become NA
  ext_tab <- table(ekey)
  amb_ext <- names(ext_tab)[ext_tab > 1]
  ext_first_idx <- match(ekey, ekey)  # first occurrence positions
  ext_map <- ext_first_idx
  names(ext_map) <- ekey

  m <- match(bkey, ekey)
  # drop ambiguous
  m[!is.na(m) & ekey[m] %in% amb_ext] <- NA_integer_

  n_matched   <- sum(!is.na(m))
  n_ambiguous <- sum(!is.na(match(bkey, amb_ext)))
  n_unmatched <- length(base_cells) - n_matched

  if (isTRUE(debug)) {
    message(
      "[scQCenrich-debug] external match: ",
      n_matched,
      " matched / ",
      length(base_cells),
      " total (",
      n_ambiguous,
      " ambiguous, ",
      n_unmatched,
      " unmatched)."
    )
  }
  list(
    map = m,
    n_matched = n_matched,
    n_ambiguous = n_ambiguous,
    n_unmatched = n_unmatched
  )
}



.dbg_file <- function()
  getOption("scQCenrich.debug_file", "qc_outputs/sctype_debug.log")
.dbg <- function(fmt, ...) {
  if (!isTRUE(getOption("scQCenrich.debug", FALSE)))
    return(invisible())
  dir.create(dirname(.dbg_file()),
             showWarnings = FALSE,
             recursive = TRUE)
  msg <- sprintf("[%s] %s\n", format(Sys.time(), "%F %T"), sprintf(fmt, ...))
  cat(msg, file = .dbg_file(), append = TRUE)
}
.capture <- function(expr, tag = "") {
  withCallingHandlers(
    expr,
    warning = function(w) {
      .dbg("[%s][WARN] %s", tag, conditionMessage(w))
      invokeRestart("muffleWarning")
    },
    message = function(m) {
      .dbg("[%s][MSG] %s", tag, conditionMessage(m))
    }
  )
}

.capture_logs <- function(expr, tag = "ScType") {
  withCallingHandlers(
    tryCatch(
      expr,
      error = function(e) {
        .dbg("%s ERROR: %s", tag, conditionMessage(e))
        NULL
      }
    ),
    warning = function(w) {
      .dbg("%s WARNING: %s", tag, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
}

#' Safe wrapper around AddModuleScore that returns the appended column names
#' @keywords internal
.add_module_score_safe <- function(obj,
                                   features,
                                   assay,
                                   name = "MS__",
                                   nbin = 24,
                                   ctrl = 100,
                                   debug_dir = NULL) {
  stopifnot(inherits(obj, "Seurat"))
  old <- colnames(obj@meta.data)
  Seurat::DefaultAssay(obj) <- assay
  obj <- Seurat::AddModuleScore(
    object = obj,
    features = features,
    assay = assay,
    nbin = nbin,
    ctrl = ctrl,
    name = name
  )
  new <- setdiff(colnames(obj@meta.data), old)

  # Order the new columns by the numeric suffix Seurat appends, e.g. MS__1, MS__2, ...
  ord <- order(suppressWarnings(as.integer(gsub(
    "^.*?([0-9]+)$", "\\1", new
  ))))
  new <- new[ord]

  # Rename to MS__<signature_name> in the SAME order as 'features'
  clean <- make.names(names(features), unique = TRUE)
  stopifnot(length(new) == length(clean))
  colnames(obj@meta.data)[match(new, colnames(obj@meta.data))] <- paste0(sub("_*$", "_", name), clean)

  if (!is.null(debug_dir)) {
    dir.create(debug_dir, FALSE, TRUE)
    write.csv(
      data.frame(
        new_col = paste0(sub("_*$", "_", name), clean),
        sig     = names(features),
        n_genes = vapply(features, length, 1L),
        stringsAsFactors = FALSE,
        check.names = FALSE
      ),
      file.path(debug_dir, "01_modulescore_columns.csv"),
      row.names = FALSE
    )
  }
  list(obj = obj, cols = paste0(sub("_*$", "_", name), clean))
}



.trace_start <- function(stage, root = file.path("qc_outputs", paste0("debug_", format(Sys.time(), "%Y%m%d-%H%M%S")))) {
  dir.create(root, FALSE, TRUE)
  cat(
    sprintf("[%s] START %s\n", Sys.time(), stage),
    file = file.path(root, "trace.log"),
    append = TRUE
  )
  list(dir = root, stage = stage)
}
.trace_log <- function(tr, fmt, ...) {
  msg <- sprintf(fmt, ...)
  message(msg)
  cat(
    sprintf("[%s] %s\n", Sys.time(), msg),
    file = file.path(tr$dir, "trace.log"),
    append = TRUE
  )
}
.trace_save <- function(tr, x, name) {
  p <- file.path(tr$dir, paste0(name, if (grepl("\\.rds$", name, TRUE))
    ""
    else
      ".rds"))
  saveRDS(x, p)
  invisible(p)
}
.trace_end <- function(tr) {
  cat(
    sprintf("[%s] END %s\n", Sys.time(), tr$stage),
    file = file.path(tr$dir, "trace.log"),
    append = TRUE
  )
}

.say <- function(...)
  if (isTRUE(getOption("scQCenrich.debug")))
    message(sprintf(...))

.top2 <- function(v) {
  o <- sort(v, decreasing = TRUE)
  c(top = o[1], second = ifelse(length(o) >= 2, o[2], -Inf))
}

.scq_dbg <- function(...) if (isTRUE(getOption("scQCenrich.debug", FALSE))) message(sprintf(...))

scq_finite <- function(x) { x[is.finite(x)] }

scq_safe_range <- function(x, pad = 0.04, fallback = c(0, 1)) {
  x <- scq_finite(x)
  if (!length(x)) return(fallback)
  r <- range(x)
  if (!all(is.finite(r))) return(fallback)
  if (r[1] == r[2]) {
    w <- if (r[1] == 0) 1 else abs(r[1]) * max(pad, 1e-6)
    r <- r + c(-w, w)
  }
  r
}

# Validate a PNG; return TRUE/FALSE
scq_is_valid_png <- function(path) {
  if (!file.exists(path)) return(FALSE)
  tryCatch({ png::readPNG(path); TRUE }, error = function(e) FALSE)
}

# Make a placeholder PNG explaining what failed
scq_write_placeholder <- function(path, title = "Plot unavailable", subtitle = NULL) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  grDevices::png(filename = path, width = 1200, height = 400, res = 150)
  old <- par(mar = c(0, 0, 2, 0)); on.exit(par(old), add = TRUE)
  plot.new(); title(main = title, cex.main = 1.1)
  if (!is.null(subtitle)) mtext(subtitle, side = 3, line = -1, cex = 0.9)
  grDevices::dev.off()
}

# Run plotting code safely; always leave behind a valid PNG
scq_png <- function(path, width = 1800, height = 1200, res = 200, expr) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  ok <- TRUE; emsg <- NULL
  grDevices::png(filename = path, width = width, height = height, res = res)
  tryCatch(
    eval.parent(substitute(expr)),
    error = function(e) { ok <<- FALSE; emsg <<- conditionMessage(e) }
  )
  grDevices::dev.off()

  if (!ok || !scq_is_valid_png(path)) {
    .scq_dbg("[scq_png] invalid PNG written to %s; err=%s", path, if (is.null(emsg)) "<none>" else emsg)
    scq_write_placeholder(path, title = "Plot failed", subtitle = emsg)
  }
}

# Debug helpers for vectors used as axes
scq_dbg_vec <- function(x, label) {
  xf <- scq_finite(x)
  .scq_dbg("[%s] N=%d finite=%d range=(%.4f, %.4f) unique=%d",
           label, length(x), length(xf),
           if (length(xf)) min(xf) else NaN,
           if (length(xf)) max(xf) else NaN,
           length(unique(xf)))
}
scq_inspect_names <- function(x, label="") {
  .scq_dbg("[names:%s] class=%s typeof=%s length=%d anyNA=%s first=%s",
           label, paste(class(x), collapse=","), typeof(x), length(x),
           any(is.na(x)), paste(utils::head(as.character(x), 3L), collapse=", "))
}

scq_make_unique <- function(x, label="") {
  if (is.null(x)) return(NULL)
  scq_inspect_names(x, paste0(label, ":raw"))
  if (!is.character(x)) {
    .scq_dbg("[names:%s] coercing to character from class=%s typeof=%s",
             label, paste(class(x), collapse=","), typeof(x))
    x <- as.character(x)
  }
  x[is.na(x)] <- "NA"
  x[x == ""]  <- "unnamed"
  # defensively drop attributes that sometimes sneak in (e.g., AsIs/levels)
  attributes(x) <- attributes(x)[intersect(names(attributes(x)), c("names"))]
  x <- base::make.unique(x)
  scq_inspect_names(x, paste0(label, ":unique"))
  x
}
