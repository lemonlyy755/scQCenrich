# ---- helpers for debug printing ----
.scqce_dbg <- function(...) {
  if (isTRUE(getOption("scQCenrich.debug", FALSE))) {
    cat(sprintf("[scQCenrich-debug][QCmetrics] %s\n", sprintf(...)))
  }
}

#' Normalize cell barcodes in several safe modes (used for cross-obj alignment)
#' @param cells character vector of cell names
#' @param mode  "basic" strips sample prefix (before ':'), strips trailing 'x',
#'              and strips trailing '-<digits>' (e.g., '-1'); "core16" returns
#'              just the last 16 A/C/G/T characters if present.
#' @keywords internal
.norm_cells <- function(cells, mode = c("basic", "core16")) {
  mode <- match.arg(mode)
  z <- cells
  if (mode == "basic") {
    # drop sample prefix "sample:BARCODE"
    z <- sub("^.*?:", "", z, perl = TRUE)
    # drop trailing x (velocyto/loom convention seen in practice)
    z <- sub("x$", "", z, perl = TRUE)
    # drop trailing -<digits> (Cell Ranger GEM group e.g. -1, -2, ...)
    z <- sub("-[0-9]+$", "", z, perl = TRUE)
    return(z)
  } else {
    # keep only a canonical 16-mer barcode core if present (A/C/G/T only)
    # take the LAST 16-mer to be permissive to prefixes/suffixes
    core <- sub(".*?([ACGT]{16}).*$", "\\1", toupper(z), perl = TRUE)
    # if no 16-mer, keep original to avoid over-normalizing
    bad <- !grepl("^[ACGT]{16}$", core)
    core[bad] <- z[bad]
    return(core)
  }
}

#' Try to align matrix columns to target cells with multiple guarded passes
#' Returns a matrix subset/reordered to target_cells; attributes include
#' alignment stats and the pass used. Never introduces regressions: if a pass
#' produces duplicate/ambiguous matches, it is skipped.
#' @keywords internal
.align_cols_robust <- function(mat, target_cells, tag = "ext") {
  stopifnot(!is.null(colnames(mat)))
  src <- colnames(mat)
  tgt <- target_cells

  pass_try <- function(src_key, tgt_key, pass_name) {
    # Build map from src_key -> src_col, ensure uniqueness on both sides
    if (anyDuplicated(src_key) || anyDuplicated(tgt_key)) {
      .scqce_dbg("[%s] pass=%s skipped due to duplicates (src_dup=%s tgt_dup=%s)",
                 tag, pass_name, anyDuplicated(src_key), anyDuplicated(tgt_key))
      return(NULL)
    }
    m <- match(tgt_key, src_key)
    ok <- !is.na(m)
    n_ok <- sum(ok)
    .scqce_dbg("[%s] pass=%s matches: %d / %d", tag, pass_name, n_ok, length(tgt))
    if (!n_ok) return(NULL)
    list(idx = m, ok = ok, name = pass_name, n_ok = n_ok)
  }

  # Pass 0: exact
  res <- pass_try(src, tgt, "exact")
  if (is.null(res) || res$n_ok < length(tgt)) {
    # Pass 1: basic normalization on both sides
    res1 <- pass_try(.norm_cells(src, "basic"), .norm_cells(tgt, "basic"), "basic")
    if (!is.null(res1) && (!res$n_ok || res1$n_ok > res$n_ok)) res <- res1
  }

  if (is.null(res) || res$n_ok == 0) {
    # Pass 2 (last resort): 16-mer core, but only if it does not collide
    src_core <- .norm_cells(src, "core16")
    tgt_core <- .norm_cells(tgt, "core16")

    # guardrail: require 1:1 mapping on the subset we’ll match
    if (!anyDuplicated(src_core) && !anyDuplicated(tgt_core)) {
      res2 <- pass_try(src_core, tgt_core, "core16")
      if (!is.null(res2)) res <- res2
    } else {
      .scqce_dbg("[%s] core16 skipped due to duplicate cores (potential collisions).", tag)
    }
  }

  if (is.null(res)) {
    .scqce_dbg("[%s] NO MATCH after all passes. example target cells: %s", tag,
               paste(utils::head(tgt, 5), collapse = ", "))
    attr(mat, "align_info") <- list(pass = "none", matched = 0L, of = length(tgt))
    return(mat[, 0, drop = FALSE]) # empty, caller will handle NA fill
  }

  # Reorder subset to target order (unmatched -> drop)
  aligned <- mat[, res$idx[res$ok], drop = FALSE]
  colnames(aligned) <- tgt[res$ok] # harmonize names to target cells

  # Write a small debug file if debugging is on
  if (isTRUE(getOption("scQCenrich.debug", FALSE))) {
    dbg <- data.frame(
      target = tgt,
      src = src[res$idx],
      matched = res$ok,
      pass = res$name,
      stringsAsFactors = FALSE
    )
    fn <- file.path("qc_outputs", sprintf("splice_alignment_debug_%s.csv", tag))
    dir.create(dirname(fn), showWarnings = FALSE, recursive = TRUE)
    utils::write.csv(dbg, fn, row.names = FALSE, quote = TRUE)
    .scqce_dbg("[%s] wrote %s", tag, fn)
  }

  attr(aligned, "align_info") <- list(pass = res$name, matched = res$n_ok, of = length(tgt))
  aligned
}

# ---- PATCH inside calcQCmetrics(): replace how you align spliced/unspliced ----
# After you obtain `spliced_mat_raw` and `unspliced_mat_raw` from external objects,
# replace your alignment code with the following block.

.align_splice_mats <- function(spliced_mat_raw, unspliced_mat_raw, target_cells) {
  # Align columns to main object cell order with guarded passes (no regression)
  sp_al <- if (!is.null(spliced_mat_raw)) .align_cols_robust(spliced_mat_raw, target_cells, tag = "spliced_ext") else NULL
  un_al <- if (!is.null(unspliced_mat_raw)) .align_cols_robust(unspliced_mat_raw, target_cells, tag = "unspliced_ext") else NULL

  # For NA-safe summaries we’ll create full-length vectors and fill aligned cols
  n_cells <- length(target_cells)

  # initialize NA vectors
  nCount_spliced   <- setNames(rep(NA_real_, n_cells), target_cells)
  nFeature_spliced <- nCount_spliced
  nCount_unspliced <- nCount_spliced
  nFeature_unspliced <- nCount_spliced

  # fill if available
  if (!is.null(sp_al) && ncol(sp_al)) {
    cs <- Matrix::colSums(sp_al)
    nf <- Matrix::colSums(sp_al > 0)
    nCount_spliced[colnames(sp_al)]   <- cs
    nFeature_spliced[colnames(sp_al)] <- nf
    info <- attr(sp_al, "align_info"); .scqce_dbg("[spliced_ext] aligned %d/%d (pass=%s)", info$matched, info$of, info$pass)
  } else {
    .scqce_dbg("[spliced_ext] aligned 0/%d", n_cells)
  }

  if (!is.null(un_al) && ncol(un_al)) {
    cs <- Matrix::colSums(un_al)
    nf <- Matrix::colSums(un_al > 0)
    nCount_unspliced[colnames(un_al)]   <- cs
    nFeature_unspliced[colnames(un_al)] <- nf
    info <- attr(un_al, "align_info"); .scqce_dbg("[unspliced_ext] aligned %d/%d (pass=%s)", info$matched, info$of, info$pass)
  } else {
    .scqce_dbg("[unspliced_ext] aligned 0/%d", n_cells)
  }

  # fractions
  tot <- nCount_spliced + nCount_unspliced
  spliced_frac   <- ifelse(tot > 0, nCount_spliced / tot, NA_real_)
  unspliced_frac <- ifelse(tot > 0, nCount_unspliced / tot, NA_real_)
  u2s_ratio      <- ifelse(nCount_spliced > 0, nCount_unspliced / nCount_spliced, NA_real_)
  intronic_frac  <- ifelse(tot > 0, nCount_unspliced / tot, NA_real_)

  list(
    nCount_spliced = nCount_spliced,
    nFeature_spliced = nFeature_spliced,
    nCount_unspliced = nCount_unspliced,
    nFeature_unspliced = nFeature_unspliced,
    spliced_frac = spliced_frac,
    unspliced_frac = unspliced_frac,
    u2s_ratio = u2s_ratio,
    intronic_frac = intronic_frac
  )
}


#' Calculate QC metrics (UMIs, genes, mito%, MALAT1, stress; optional splicing metrics)
#' Computes nCount/nFeature, pctMT, MALAT1 fraction, stress_score (case-insensitive SYMBOL match),
#' and optional splicing metrics from same object layers or external Seurat objects.
#' Compatible with Seurat v5 layers and v4 slots. Returns a data.frame and (if add_to_meta=TRUE) writes
#' *_QC fields into obj@meta.data (no regression in field names).
#' @param obj Seurat object
#' @param assay assay name (default "RNA")
#' @param slot for v4 compatibility ("counts"); ignored when layer is used internally
#' @param species "mouse" or "human"
#' @param spliced,unspliced optional layer spec within same object; either a single layer string (e.g., "spliced")
#'        or "assay:layer" (e.g., "spliced:counts"). Ignored if external objects are provided.
#' @param spliced_obj,unspliced_obj optional external Seurat objects
#' @param spliced_assay,unspliced_assay assay names on external objects
#' @param spliced_layer,unspliced_layer layer names on external objects (e.g., "counts")
#' @param add_to_meta write *_QC back to obj@meta.data
#' @param debug logical; defaults to options(scQCenrich.debug, FALSE)
#' @export
calcQCmetrics <- function(
    obj,
    assay = "RNA",
    slot  = "counts",              # v4 fallback only
    species = c("mouse","human"),
    spliced = NULL, unspliced = NULL,             # same-object layers, or:
    spliced_obj = NULL, spliced_assay = NULL, spliced_layer = "counts",
    unspliced_obj = NULL, unspliced_assay = NULL, unspliced_layer = "counts",
    add_to_meta = TRUE,
    debug = getOption("scQCenrich.debug", FALSE)
){
  stopifnot(inherits(obj, "Seurat"))
  species <- match.arg(species)
  cells <- colnames(obj)

  # ---- tiny logging helper ----
  .dbg <- function(...){
    if (isTRUE(debug)) message("[scQCenrich-debug][QCmetrics] ", sprintf(...))
  }

  # ---- helpers ----
  .as_dgc <- function(x, nm="matrix"){
    if (inherits(x, "dgCMatrix")) return(x)
    if (is.null(x)) stop(sprintf("[calcQCmetrics] %s is NULL", nm))
    if (is.matrix(x)) return(Matrix::Matrix(x, sparse = TRUE))
    y <- tryCatch(as.matrix(x), error = function(e) NULL)
    if (!is.null(y)) return(Matrix::Matrix(y, sparse = TRUE))
    stop(sprintf("[calcQCmetrics] could not coerce %s to dgCMatrix", nm))
  }

  .get_layer <- function(object, assay, layer, fallback_slot=NULL){
    # Seurat v5 prefers layers; slot is fallback (v4) :contentReference[oaicite:1]{index=1}
    mat <- tryCatch(SeuratObject::GetAssayData(object=object, assay=assay, layer=layer),
                    error=function(e) NULL)
    if (!is.null(mat) && length(dim(mat))==2L) {
      .dbg("GetAssayData(layer='%s') -> %s x %s", layer, nrow(mat), ncol(mat))
      return(.as_dgc(mat, sprintf("%s:%s [layer]", assay, layer)))
    }
    if (!is.null(fallback_slot)) {
      mat2 <- tryCatch(SeuratObject::GetAssayData(object=object, assay=assay, slot=fallback_slot),
                       error=function(e) NULL)
      if (!is.null(mat2) && length(dim(mat2))==2L) {
        .dbg("GetAssayData(slot='%s') -> %s x %s", fallback_slot, nrow(mat2), ncol(mat2))
        return(.as_dgc(mat2, sprintf("%s:%s [slot]", assay, fallback_slot)))
      }
    }
    aa <- tryCatch(object[[assay]], error=function(e) NULL)
    if (!is.null(aa) && !is.null(aa@layers) && !is.null(aa@layers[[layer]])) {
      mat3 <- aa@layers[[layer]]
      .dbg("assay@layers[['%s']] -> %s x %s", layer, nrow(mat3), ncol(mat3))
      return(.as_dgc(mat3, sprintf("%s:%s [layers]", assay, layer)))
    }
    .dbg("Failed to obtain %s:%s (if multiple layers exist, consider JoinLayers()).", assay, layer) # :contentReference[oaicite:2]{index=2}
    stop(sprintf("[calcQCmetrics] Failed to obtain %s:%s (layer/slot).", assay, layer))
  }

  .safe_colsums <- function(m){
    if (is.null(m) || length(m)==0L) return(numeric(0))
    if (is.null(dim(m))) m <- Matrix::Matrix(matrix(m, nrow=1L), sparse=TRUE)
    if (nrow(m)==0L) return(rep(0, ncol(m)))
    Matrix::colSums(m)
  }

  .match_layer_same_obj <- function(spec){
    if (is.null(spec)) return(NULL)
    if (grepl(":", spec, fixed=TRUE)) {
      pr <- strsplit(spec, ":", fixed=TRUE)[[1]]
      .get_layer(obj, pr[1], pr[2], fallback_slot = pr[2])
    } else {
      .get_layer(obj, assay, spec, fallback_slot = spec)
    }
  }

  # ---- barcode normalizer (your variant) + guarded 16-mer fallback ----
  .norm_cells <- function(v, mode=c("basic","core16")){
    mode <- match.arg(mode)
    z <- v
    if (mode == "basic") {
      # drop sample prefix before ':'; drop trailing 'x'; drop trailing '-<digits>'
      z <- sub("^.*?:", "", z, perl=TRUE)
      z <- sub("x$", "", z, perl=TRUE)
      z <- sub("-[0-9]+$", "", z, perl=TRUE)
      return(z)
    } else {
      core <- sub(".*?([ACGT]{16}).*$", "\\1", toupper(z), perl=TRUE)  # 10x barcodes are 16 bp :contentReference[oaicite:3]{index=3}
      bad <- !grepl("^[ACGT]{16}$", core)
      core[bad] <- z[bad]
      return(core)
    }
  }

  .align_cols_robust <- function(mat, target_cells, tag="ext"){
    stopifnot(!is.null(colnames(mat)))
    src <- colnames(mat); tgt <- target_cells

    pass_try <- function(src_key, tgt_key, pass_name){
      if (anyDuplicated(src_key) || anyDuplicated(tgt_key)) {
        .dbg("[%s] pass=%s skipped (duplicates in keys)", tag, pass_name)
        return(NULL)
      }
      idx <- match(tgt_key, src_key)
      ok  <- !is.na(idx); n_ok <- sum(ok)
      .dbg("[%s] pass=%s matches: %d / %d", tag, pass_name, n_ok, length(tgt))
      if (!n_ok) return(NULL)
      list(idx=idx, ok=ok, name=pass_name, n_ok=n_ok)
    }

    # Pass 0: exact
    best <- pass_try(src, tgt, "exact")

    # Pass 1: your basic normalizer on both sides
    alt1 <- pass_try(.norm_cells(src,"basic"), .norm_cells(tgt,"basic"), "basic")
    if (is.null(best) || (!is.null(alt1) && alt1$n_ok > best$n_ok)) best <- alt1

    # Pass 2: core16 (only if unique)
    src16 <- .norm_cells(src,"core16"); tgt16 <- .norm_cells(tgt,"core16")
    if (!anyDuplicated(src16) && !anyDuplicated(tgt16)) {
      alt2 <- pass_try(src16, tgt16, "core16")
      if (is.null(best) || (!is.null(alt2) && alt2$n_ok > best$n_ok)) best <- alt2
    } else {
      .dbg("[%s] core16 skipped due to duplicate cores (potential collisions).", tag)
    }

    if (is.null(best)) {
      .dbg("[%s] NO MATCH after all passes. example targets: %s", tag, paste(utils::head(tgt,5), collapse=", "))
      attr(mat, "align_info") <- list(pass="none", matched=0L, of=length(tgt))
      return(mat[, 0, drop=FALSE]) # empty -> caller fills NA
    }

    aligned <- mat[, best$idx[best$ok], drop=FALSE]
    colnames(aligned) <- tgt[best$ok]
    if (isTRUE(debug)) {
      dir.create("qc_outputs", showWarnings=FALSE, recursive=TRUE)
      dbg <- data.frame(target=tgt, src=src[best$idx], matched=best$ok, pass=best$name, stringsAsFactors=FALSE)
      utils::write.csv(dbg, file=file.path("qc_outputs", sprintf("splice_alignment_debug_%s.csv", tag)), row.names=FALSE)
      .dbg("[%s] wrote qc_outputs/splice_alignment_debug_%s.csv", tag, tag)
    }
    attr(aligned, "align_info") <- list(pass=best$name, matched=best$n_ok, of=length(tgt))
    aligned
  }

  # ---- base counts & basic QC ----
  counts <- .get_layer(obj, assay, layer="counts", fallback_slot="counts")
  counts <- counts[, cells, drop=FALSE]
  .dbg("Base counts dim: %s x %s (assay=%s)", nrow(counts), ncol(counts), assay)

  nCount   <- .safe_colsums(counts)
  nFeature <- Matrix::colSums(counts > 0)

  rn <- rownames(counts)
  mito_pat <- if (species=="mouse") "^mt-" else "^MT-"
  pctMT <- if (any(grepl(mito_pat, rn, ignore.case=TRUE))) {
    Matrix::colSums(counts[grepl(mito_pat, rn, ignore.case=TRUE), , drop=FALSE]) / pmax(nCount,1)
  } else rep(NA_real_, length(nCount))
  .dbg("Mito pattern='%s' | mito genes found: %s", mito_pat, sum(grepl(mito_pat, rn, ignore.case=TRUE)))

  malat_idx <- which(toupper(rn) %in% toupper(c("MALAT1","Malat1")))
  MALAT1_frac <- if (length(malat_idx)) Matrix::colSums(counts[malat_idx, , drop=FALSE]) / pmax(nCount,1) else rep(NA_real_, length(nCount))
  .dbg("MALAT1 rows matched: %s", length(malat_idx))

  # normalized (for stress score)
  norm <- tryCatch(.get_layer(obj, assay, layer="data", fallback_slot="data"), error=function(e) NULL)
  if (is.null(norm) || nrow(norm)==0L){
    lib <- pmax(nCount,1)
    norm <- log1p(t(t(counts)/lib) * 1e6)
    norm <- .as_dgc(norm, "logCPM")
    .dbg("Norm layer missing -> using logCPM fallback (rows=%s, cols=%s)", nrow(norm), ncol(norm))
  } else {
    norm <- norm[, cells, drop=FALSE]
    .dbg("Norm data (layer='data') dims: %s x %s", nrow(norm), ncol(norm))
  }

  # ---- spliced/unspliced (same object or external) ----
  get_ext <- function(o, a, l, label){
    .dbg("[%s] trying assay=%s layer/slot=%s", label, a, l)
    m <- .get_layer(o, a, layer=l, fallback_slot=l)
    .dbg("[%s] raw dim: %s x %s", label, nrow(m), ncol(m))
    m
  }

  spliced_raw <- if (!is.null(spliced_obj)) {
    stopifnot(inherits(spliced_obj,"Seurat"), !is.null(spliced_assay))
    get_ext(spliced_obj, spliced_assay, spliced_layer, "spliced_ext")
  } else if (!is.null(spliced)) .match_layer_same_obj(spliced) else NULL

  unspliced_raw <- if (!is.null(unspliced_obj)) {
    stopifnot(inherits(unspliced_obj,"Seurat"), !is.null(unspliced_assay))
    get_ext(unspliced_obj, unspliced_assay, unspliced_layer, "unspliced_ext")
  } else if (!is.null(unspliced)) .match_layer_same_obj(unspliced) else NULL

  # Align columns to the main object’s cell order (safe passes only)
  sp_al <- if (!is.null(spliced_raw)) .align_cols_robust(spliced_raw, cells, tag="spliced_ext") else NULL
  un_al <- if (!is.null(unspliced_raw)) .align_cols_robust(unspliced_raw, cells, tag="unspliced_ext") else NULL

  .dbg("Aligned spliced cols: %s / %s | unspliced cols: %s / %s",
       if (!is.null(sp_al)) ncol(sp_al) else 0L, length(cells),
       if (!is.null(un_al)) ncol(un_al) else 0L, length(cells))

  # ---- compute splicing metrics (fill NA when not aligned) ----
  nCount_spliced   <- setNames(rep(NA_real_, length(cells)), cells)
  nFeature_spliced <- nCount_spliced
  nCount_unspliced <- nCount_spliced
  nFeature_unspliced <- nCount_spliced

  if (!is.null(sp_al) && ncol(sp_al)>0) {
    nCount_spliced[colnames(sp_al)]   <- Matrix::colSums(sp_al)
    nFeature_spliced[colnames(sp_al)] <- Matrix::colSums(sp_al > 0)
    info <- attr(sp_al, "align_info"); .dbg("[spliced_ext] aligned %d/%d (pass=%s)", info$matched, info$of, info$pass)
  } else .dbg("[spliced_ext] aligned 0/%d", length(cells))

  if (!is.null(un_al) && ncol(un_al)>0) {
    nCount_unspliced[colnames(un_al)]   <- Matrix::colSums(un_al)
    nFeature_unspliced[colnames(un_al)] <- Matrix::colSums(un_al > 0)
    info <- attr(un_al, "align_info"); .dbg("[unspliced_ext] aligned %d/%d (pass=%s)", info$matched, info$of, info$pass)
  } else .dbg("[unspliced_ext] aligned 0/%d", length(cells))

  tot <- nCount_spliced + nCount_unspliced
  spliced_frac   <- ifelse(tot > 0, nCount_spliced / tot, NA_real_)
  unspliced_frac <- ifelse(tot > 0, nCount_unspliced / tot, NA_real_)
  u2s_ratio      <- ifelse(nCount_spliced > 0, nCount_unspliced / nCount_spliced, NA_real_)
  intronic_frac  <- ifelse(tot > 0, nCount_unspliced / tot, NA_real_)
  .dbg("Cells with tot(spliced+unspliced) <= 0: %s (these yield NA intronic_frac).", sum(!(tot > 0), na.rm=TRUE))

  # ---- stress score (SYMBOL, case-insensitive) ----
  stress_genes <- default_stress_genes(species)  # SYMBOLs
  sg <- rownames(norm)[toupper(rownames(norm)) %in% toupper(stress_genes)]
  stress_score <- if (length(sg)) Matrix::colMeans(norm[sg, , drop=FALSE]) else rep(NA_real_, ncol(norm))
  .dbg("Stress genes matched in norm matrix: %s", length(sg))

  # ---- assemble & (optionally) write back ----
  out <- data.frame(
    nCount = as.numeric(nCount),
    nFeature = as.numeric(nFeature),
    pctMT = as.numeric(pctMT),
    MALAT1_frac = as.numeric(MALAT1_frac),
    stress_score = as.numeric(stress_score),
    intronic_frac = as.numeric(intronic_frac),
    nCount_spliced = as.numeric(nCount_spliced[cells]),
    nFeature_spliced = as.numeric(nFeature_spliced[cells]),
    nCount_unspliced = as.numeric(nCount_unspliced[cells]),
    nFeature_unspliced = as.numeric(nFeature_unspliced[cells]),
    spliced_frac = as.numeric(spliced_frac),
    unspliced_frac = as.numeric(unspliced_frac),
    u2s_ratio = as.numeric(u2s_ratio),
    row.names = cells,
    check.names = FALSE
  )

  # keep doublet meta if present
  md <- obj@meta.data
  out$is_doublet <- if ("is_doublet" %in% colnames(md)) as.logical(md[cells,"is_doublet"]) else FALSE
  out$dbl_score  <- if ("dbl_score"  %in% colnames(md)) as.numeric(md[cells,"dbl_score"])  else NA_real_

  if (isTRUE(add_to_meta)) {
    md[cells, "nCount_QC"]         <- out$nCount
    md[cells, "nFeature_QC"]       <- out$nFeature
    md[cells, "pctMT_QC"]          <- out$pctMT
    md[cells, "MALAT1_frac_QC"]    <- out$MALAT1_frac
    md[cells, "stress_score_QC"]   <- out$stress_score
    md[cells, "intronic_frac_QC"]  <- out$intronic_frac
    md[cells, "spliced_frac_QC"]   <- out$spliced_frac
    md[cells, "unspliced_frac_QC"] <- out$unspliced_frac
    md[cells, "u2s_ratio_QC"]      <- out$u2s_ratio
    obj@meta.data <- md
  }

  attr(out, "obj") <- obj
  out
}





#' Default stress genes
#' @param species 'mouse' or 'human'
#' @return Character vector of stress-related genes
#' @export
default_stress_genes <- function(species = c("mouse","human")){
  species <- match.arg(species)
  if (species == "mouse") {
    unique(c("Fos","Fosb","Jun","Junb","Jund","Atf3","Egr1","Dusp1","Dusp2",
             "Ier2","Ier3","Btg1","Btg2","Zfp36","Zfp36l1","Hspa8","Hsp90ab1",
             "Hspa1a","Hspa1b","Hspb1","Hsp90aa1","Hsp90b1","Hsph1","Hspe1",
             "Hspd1","Xbp1","Ddit3"))
  } else {
    unique(c("FOS","FOSB","JUN","JUNB","JUND","ATF3","EGR1","DUSP1","DUSP2",
             "IER2","IER3","BTG1","BTG2","ZFP36","ZFP36L1","HSPA8","HSP90AB1",
             "HSPA1A","HSPA1B","HSPB1","HSP90AA1","HSP90B1","HSPH1","HSPE1",
             "HSPD1","XBP1","DDIT3"))
  }
}
