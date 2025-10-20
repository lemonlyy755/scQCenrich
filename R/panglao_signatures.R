
panglao_signatures <- function(
    tsv,
    species = c("human","mouse"),
    canonical_only = TRUE,
    ui_max = 0.20,
    tissue = NULL,
    min_genes = 5,
    debug = getOption("scQCenrich.debug", FALSE),
    debug_dir = file.path("qc_outputs","panglao_debug"),
    obj = NULL,
    assay = if (inherits(obj,"Seurat")) Seurat::DefaultAssay(obj) else NULL
){
  species <- match.arg(species)
  if (!requireNamespace("readr", quietly = TRUE)) stop("Need 'readr' to read Panglao TSV")

  .say <- function(...) if (debug) message(sprintf(...))
  stopifnot(file.exists(tsv))
  df <- readr::read_tsv(tsv, show_col_types = FALSE)

  # --- normalize headers & map likely column names
  cn0 <- names(df); cn <- tolower(gsub("\\s+","_", cn0)); names(df) <- cn
  pick <- function(cands) { x <- intersect(tolower(cands), names(df)); if (length(x)) x[1] else NA_character_ }
  col_gene    <- pick(c("official_gene_symbol","gene_symbol","gene"))
  col_ct      <- pick(c("cell_type","celltype","cell"))
  col_organ   <- pick(c("organ","tissue","tissue_type","tissue_name"))  # <-- organ/tissue column
  col_species <- pick(c("species","organism","org"))
  col_canon   <- pick(c("canonical_marker","canonical"))
  col_ui      <- pick(c("ubiquitousness_index","ui"))

  if (any(is.na(c(col_gene, col_ct)))) stop("Panglao TSV must have gene symbol & cell-type columns.")

  # --- species filter (Panglao uses 'Hs'/'Mm' in many dumps)
  if (!is.na(col_species) && col_species %in% names(df)) {
    spv <- tolower(as.character(df[[col_species]]))
    is_mouse <- grepl("mm", spv, ignore.case = TRUE)
    is_human <- grepl("hs", spv, ignore.case = TRUE)
    keep <- if (species == "mouse") is_mouse | (!any(is_mouse, na.rm=TRUE) & !any(is_human, na.rm=TRUE))
    else                     is_human | (!any(is_mouse, na.rm=TRUE) & !any(is_human, na.rm=TRUE))
    df <- df[ifelse(is.finite(keep), keep, TRUE), , drop = FALSE]
  }

  # --- canonical-only filter (if available)
  if (isTRUE(canonical_only) && !is.na(col_canon) && col_canon %in% names(df)) {
    canon_raw <- df[[col_canon]]
    is_canon  <- canon_raw %in% c(1L,"1",TRUE,"TRUE","True")
    df <- df[is_canon %in% TRUE, , drop = FALSE]
  }

  # --- UI filter (specificity), if present
  if (!is.null(ui_max) && !is.na(col_ui) && col_ui %in% names(df)) {
    ui <- suppressWarnings(as.numeric(df[[col_ui]]))
    df <- df[ is.na(ui) | ui <= ui_max, , drop = FALSE]
  }

  # --- NEW: tissue/organ filter (vector-friendly; case-insensitive; partial match fallback)
  if (!is.null(tissue) && length(tissue) && !is.na(col_organ) && col_organ %in% names(df)) {
    orgv <- as.character(df[[col_organ]])
    n_before <- nrow(df)
    q <- tolower(trimws(as.character(tissue)))
    keep <- tolower(trimws(orgv)) %in% q
    if (!any(keep, na.rm = TRUE)) {
      pat <- paste0(gsub("\\s+","\\\\s*", q), collapse = "|")  # partial match across any provided tissues
      keep <- grepl(pat, tolower(orgv))
    }
    if (any(keep, na.rm = TRUE)) {
      df <- df[keep %in% TRUE, , drop = FALSE]
      .say("[Panglao] tissue filter kept %d/%d rows", nrow(df), n_before)
    } else {
      .say("[Panglao] tissue filter matched 0 rows; skipping tissue filter.")
    }
  }

  # --- build signatures (uppercase symbols) & tidy
  genes <- toupper(as.character(df[[col_gene]]))
  ctype <- as.character(df[[col_ct]])
  sig   <- split(genes, ctype)
  sig   <- lapply(sig, function(v) unique(v[nzchar(v)]))

  # keep sufficient-sized sets
  sig <- sig[vapply(sig, length, integer(1)) >= min_genes]
  sig
}



panglao_signatures_debug <- function(
    tsv, species = c("human","mouse"),
    canonical_only = TRUE, ui_max = 0.20, tissue = NULL, min_genes = 5,
    obj = NULL, assay = if (inherits(obj,"Seurat")) Seurat::DefaultAssay(obj) else NULL,
    debug_dir = file.path("qc_outputs","panglao_debug")
){
  tr <- .trace_start("panglao", debug_dir)
  on.exit(.trace_end(tr), add = TRUE)
  species <- match.arg(species)

  sig <- panglao_signatures(
    tsv = tsv, species = species, canonical_only = canonical_only,
    ui_max = ui_max, tissue = tissue, min_genes = min_genes,
    debug = TRUE, debug_dir = tr$dir, obj = obj, assay = assay
  )

  .trace_log(tr, "[Panglao] signatures kept: %d (min_genes=%d, canonical_only=%s, ui_max=%s)",
             length(sig), min_genes, canonical_only, as.character(ui_max))
  .trace_save(tr, sig, "panglao_signatures.rds")

  # coverage vs object
  if (inherits(obj, "Seurat")) {
    genes_obj <- rownames(Seurat::GetAssayData(obj, assay = assay, layer = "data"))
    if (!length(genes_obj)) genes_obj <- rownames(Seurat::GetAssayData(obj, assay = assay, layer = "counts"))
    up <- toupper(genes_obj)
    tab <- data.frame(
      set      = names(sig),
      n_markers= vapply(sig, length, 1L),
      overlap  = vapply(sig, function(v) sum(toupper(v) %in% up), 1L),
      stringsAsFactors = FALSE, check.names = FALSE
    )
    tab$frac <- with(tab, ifelse(n_markers>0, overlap/n_markers, 0))
    write.csv(tab[order(-tab$frac,-tab$overlap),], file.path(tr$dir,"06_overlap_vs_object.csv"), row.names = FALSE)
    .trace_log(tr, "[Panglao] median overlap=%.2f", stats::median(tab$frac, na.rm = TRUE))
  }
  sig
}
