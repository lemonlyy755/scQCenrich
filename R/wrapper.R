




#' High-level wrapper that computes QC metrics, flags low-quality cells,
#' (optionally) detects doublets and performs lightweight annotation, and
#' renders an HTML report.
#'
#' @param obj Seurat object
#' @param species 'mouse' or 'human'
#' @param method 'threshold', 'gmm' or 'auto'
#' @param score_cutoff Integer cutoff used by 'threshold' mode
#' @param auto_k Ignored; kept for back-compat
#' @param assay Assay name
#' @param save_dir Directory for intermediate outputs
#' @param report_file Output HTML report path
#' @param spliced,unspliced Back-compat layer strings on the same object
#' @param spliced_obj,unspliced_obj Optional external Seurat objects
#' @param spliced_assay,unspliced_assay Assay names on external objects
#' @param spliced_layer,unspliced_layer Layer/slot on external objects
#' @param doublets 'auto','doubletfinder','scdblfinder','none'
#' @param remove_doublets Logical; if TRUE, mark doublets to remove/borderline
#' @param sample_col Optional meta column with sample/group
#' @param annot_method 'SingleR','marker_score','none'
#' @param marker_source 'panglao','internal','custom' (only used when annot_method='marker_score')
#' @param tissue Character vector of target organs/tissues, e.g. c("Liver","Kidney").
#' @param panglao_file Optional path to Panglao TSV (defaults to inst/extdata if present)
#' @param canonical_only,ui_max,min_genes Panglao signature filtering knobs
#' @param report_html Logical; if TRUE, render the report
#' @param theme_bootswatch Bootswatch theme for report
#' @param marker_method Either "module_score" (default) or "findmarkers".
#' @param debug Logical for verbose messages
#' @param sctype_db_path Path to ScType Excel DB. Defaults to the packaged file
#' @param annot_args List of extra args forwarded to `auto_annotate()`.
#' @param unknown_min_genes Integer; guard for "unknown" assignment.
#' @param unknown_min_margin Numeric; guard for "unknown" assignment.
#' @param unknown_min_top Integer; guard for "unknown" assignment.
#' @param qc_strength One of \code{c("lenient","moderate","strict")}; passed to \code{flagLowQuality()}.
#' @param enrichment_plots logical. If TRUE, run GO/KEGG enrichment plots; if FALSE, skip. Default: TRUE.
#' @return A list with at least: \code{$obj} (QC-kept), \code{$metrics}, \code{$status_df}, \code{$report_file}.
#'         Also returns \code{$obj_all} for convenience.
#' @export
run_qc_pipeline <- function(obj,
                            species = c("mouse", "human"),
                            method  = c("gmm","threshold",  "auto"),
                            score_cutoff = 2,
                            auto_k = 3,
                            assay = "RNA",
                            save_dir = "qc_outputs",
                            report_file = "qc_report.html",
                            # layer strings still supported
                            spliced   = NULL,
                            unspliced = NULL,
                            # external Seurat(s) with spliced/unspliced counts
                            spliced_obj = NULL,
                            spliced_assay = NULL,
                            spliced_layer = "counts",
                            unspliced_obj = NULL,
                            unspliced_assay = NULL,
                            unspliced_layer = "counts",
                            doublets  = c( "none","auto", "doubletfinder", "scdblfinder"),
                            remove_doublets = TRUE,
                            sample_col = NULL,
                            annot_method = c("marker_score", "singleR","sctype",  "none"),
                            marker_source = c("panglao", "internal", "custom"),
                            panglao_file = NULL,
                            canonical_only = TRUE,
                            ui_max = 0.20,
                            tissue = NULL,
                            min_genes = 5,
                            report_html = TRUE,
                            theme_bootswatch = "minty",
                            debug = FALSE,
                            sctype_db_path = system.file("extdata", "ScTypeDB_full.xlsx",
                                                         package = "scQCenrich", mustWork = FALSE),
                            marker_method = c("findmarkers","module_score" ) ,
                            annot_args     = NULL,
                            unknown_min_genes = 2,
                            unknown_min_margin = 0.01,
                            unknown_min_top = 0.01,
                            qc_strength = c("auto","default","lenient","strict"),
                            enrichment_plots = T
                            ) {
  species       <- match.arg(species)
  method        <- match.arg(method)
  doublets      <- match.arg(doublets)
  marker_method <- match.arg(marker_method)
  qc_strength <- match.arg(qc_strength)
  annot_method <- annot_method[1L]
  # accept both "singleR" and "SingleR" (docs mention the latter)
  annot_method  <- match.arg(annot_method,
                             c("sctype", "marker_score", "singleR", "SingleR", "none"))
  marker_source <- match.arg(marker_source)
  am <- tolower(annot_method)


  if (identical(sctype_db_path, "") || !file.exists(sctype_db_path)) {
    if (isTRUE(debug)) message("[scQCenrich-debug] No packaged ScType DB found; using NULL.")
    sctype_db_path <- NULL
  }

  o_old <- getOption("scQCenrich.debug")
  on.exit(options(scQCenrich.debug = o_old), add = TRUE)
  options(scQCenrich.debug = isTRUE(debug))

  if (isTRUE(debug)) {
    options(
      scQCenrich.debug = TRUE,
      scQCenrich.debug_file = file.path(save_dir, "sctype_debug.log")
    )
    .dbg("Debug logging to %s", getOption("scQCenrich.debug_file"))
  }

  if (!dir.exists(save_dir))
    dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)

  # [1/7] Assay/layer prep (ensure non-empty 'data' & join v5 layers if needed)
  Seurat::DefaultAssay(obj) <- assay
  obj <- .ensure_layers_and_data(obj, assay = assay)


  # [2/7] Optional annotation (cluster-level ScType / marker_score / SingleR)
  if (!identical(annot_method, "none")) {
    message(sprintf("[2/7] Auto-annotating cell types (%s) ...", annot_method))
    if (am == "singler") {
      obj <- auto_annotate(
        obj,
        species = species,
        assay = assay,
        prefer = c("singleR", "label_transfer", "marker_score"),
        level  = "cluster",
        marker_score_level = "cluster"
      )
    } else if (am == "sctype") {
      obj <- auto_annotate(
        obj,
        species = species,
        assay = assay,
        prefer = "sctype",
        level = "cluster",
        organ = tissue,
        sctype_db_path = sctype_db_path
      )  # no restraint

    } else if (am == "marker_score") {
      # choose signature source
      sig <- NULL
      if (marker_source == "panglao") {
        if (is.null(panglao_file)) {
          panglao_file <- system.file(
            "extdata",
            "PanglaoDB_markers_27_Mar_2020.tsv",
            package = utils::packageName()
          )
        }
        if (!nzchar(panglao_file) ||
            !file.exists(panglao_file)) {
          stop("Panglao file not found. Provide a valid path via `panglao_file`.")
        }
        sig <- panglao_signatures(
          tsv = panglao_file,
          species = species,
          canonical_only = canonical_only,
          ui_max = ui_max,
          tissue = tissue,
          min_genes = min_genes
        )
      } else if (marker_source == "custom") {
        stop(
          "marker_source='custom' requires passing custom_signatures to auto_annotate()."
        )
      }
      # use auto_annotate's marker_score backend (generic if sig=NULL





      obj <- auto_annotate(
        obj,


        species = species,


        assay = assay,


        prefer = "marker_score",


        level  = "cluster",


        marker_score_level = "cluster",


        custom_signatures = sig,
        fm_padj_max = 1.1,
        overlap_mode = c( "jaccard"),
        unknown_min_genes = unknown_min_genes,
        unknown_min_margin = unknown_min_margin,
        unknown_min_top = unknown_min_top,
      )
    }
  } else {
    message("[2/7] Annotation skipped (annot_method='none').")

  }

  if(enrichment_plots){
      validation_plots2_post_annotation(
    obj,
    save_dir       = "qc_outputs",
    species        = species,
    # or "human"
    assay          = assay,
    panglao_ui_max = 0.20,
    # keep canonical, specific markers
    fm_top_n       = 50,
    fm_padj_max    = 0.05,
    per_cluster_keep = 10,
    overlap_mode   = "count"
  )


  }



  # [2/7] Assay/layer prep (ensure non-empty 'data' & v5 layers joined if needed)
  Seurat::DefaultAssay(obj) <- assay
  obj <- .ensure_layers_and_data(obj, assay = assay)

  # [3/7] Doublet detection (optional)
  if (identical(doublets, "none")) {
    message("[3/7] Skipping doublet detection (disabled).")
    obj@meta.data$is_doublet <- FALSE
    obj@meta.data$dbl_score  <- NA_real_
    obj@meta.data$dbl_method <- NA_character_
  } else {
    message(sprintf("[3/7] Detecting doublets (%s) ...", doublets))
    obj@meta.data$is_doublet <- FALSE
    obj@meta.data$dbl_score  <- NA_real_
    obj@meta.data$dbl_method <- NA_character_
    dbl <- NULL
    if (doublets %in% c("doubletfinder", "auto")) {
      dbl <- detect_doublets_doubletfinder(obj, assay = assay, debug = debug)
    }
    if ((is.null(dbl) ||
         inherits(dbl, "try-error")) &&
        doublets %in% c("scdblfinder", "auto")) {
      dbl <- detect_doublets_scdblfinder(obj, assay = assay)
    }
    if (!is.null(dbl) && !inherits(dbl, "try-error") && nrow(dbl)) {
      rn <- intersect(rownames(obj@meta.data), rownames(dbl))
      obj@meta.data[rn, "is_doublet"] <- as.logical(dbl[rn, "is_doublet"])
      if ("dbl_score"  %in% colnames(dbl))
        obj@meta.data[rn, "dbl_score"]  <- as.numeric(dbl[rn, "dbl_score"])
      if ("dbl_method" %in% colnames(dbl))
        obj@meta.data[rn, "dbl_method"] <- as.character(dbl[rn, "dbl_method"])
    } else {
      warning("Doublet detection failed; proceeding without doublet flags.")
    }
  }

  # [4/7] Metrics
  message("[4/7] Calculating QC metrics ...")
  metrics <- calcQCmetrics(
    obj,
    assay = assay,
    slot = "counts",
    species = species,
    spliced = spliced,
    unspliced = unspliced,
    spliced_obj = spliced_obj,
    spliced_assay = spliced_assay,
    spliced_layer = spliced_layer,
    unspliced_obj = unspliced_obj,
    unspliced_assay = unspliced_assay,
    unspliced_layer = unspliced_layer,
    add_to_meta = TRUE
  )
  obj <- attr(metrics, "obj")

  # [5/7] Flagging (with rescue); map remove_doublets -> doublet_action
  message(sprintf(
    "[5/7] Flagging low-quality cells (%s + rescue) ...",
    if (identical(method, "auto"))
      "gmm"
    else
      method
  ))
  method2 <- if (identical(method, "auto"))
    "gmm"
  else
    method
  flags <- flagLowQuality(
    metrics,
    method       = method2,
    score_cutoff = score_cutoff,
    obj          = obj,
    celltype_col = "auto_celltype",
    sample_col   = sample_col,
    rescue_mode  = "moderate",
    qc_strength = qc_strength,
    min_cluster_size = NULL,
    save_dir = save_dir,
    doublet_action = if (isTRUE(remove_doublets))
      "remove"
    else
      "none",
    doublet_col     = "is_doublet"
  )
  stopifnot(is.data.frame(flags))
  status <- flags[rownames(metrics), , drop = FALSE]
  if (getOption("scQCenrich.debug", FALSE)) {
    outdir <- save_dir
    if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

    # full status snapshot
    utils::write.csv(status, file.path(outdir, "qc_status.csv"), row.names = TRUE)

    # convenient barcode files
    keep_ids <- rownames(status)[status$qc_status %in% c("keep","borderline")]
    rem_ids  <- rownames(status)[status$qc_status %in% c("remove")]
    writeLines(keep_ids, file.path(outdir, "barcodes_kept.txt"))
    writeLines(rem_ids,  file.path(outdir, "barcodes_removed.txt"))

    # if rescue gate table bubbled up, persist it too
    gt <- attr(flags, "rescue_gate_table")
    if (!is.null(gt)) utils::write.csv(gt, file.path(outdir, "rescue_gate_table.csv"), row.names = FALSE)

    message(sprintf("[subset] keep=%d | borderline=%d | remove=%d",
                    length(keep_ids) - sum(status$qc_status=="borderline", na.rm=TRUE),
                    sum(status$qc_status=="borderline", na.rm=TRUE),
                    length(rem_ids)))
  }

  if (isTRUE(getOption("scQCenrich.debug", FALSE))) {
    dir.create(save_dir, showWarnings = FALSE, recursive = TRUE)
    utils::write.csv(status, file.path(save_dir, "qc_status.csv"), row.names = TRUE)

    gt <- attr(flags, "rescue_gate_table")
    if (!is.null(gt)) {
      utils::write.csv(gt, file.path(save_dir, "rescue_gate_table.csv"), row.names = FALSE)
    }
    gt_alt <- attr(flags, "rescue_gate_alt")
    if (!is.null(gt_alt)) {
      utils::write.csv(gt_alt, file.path(save_dir, "rescue_gate_table_lenient.csv"), row.names = FALSE)
    }
    if (!is.null(attr(flags, "rescue_alt"))) {
      writeLines(names(which(attr(flags, "rescue_alt")$rescue_flag)),
                 con = file.path(save_dir, "rescue_lenient_candidates.txt"))
    }

    # quick console preview
    message(sprintf("[subset] keep=%d | borderline=%d | remove=%d",
                    sum(status$qc_status=="keep", na.rm=TRUE),
                    sum(status$qc_status=="borderline", na.rm=TRUE),
                    sum(status$qc_status=="remove", na.rm=TRUE)))
  }

  # robust Seurat v5-safe filter
  obj_keep <- applyQCfilter(obj, status)
  # write meta
  obj@meta.data$qc_status     <- status$qc_status
  obj@meta.data$qc_score      <- status$qc_score
  obj@meta.data$qc_cluster    <- status$cluster
  obj@meta.data$qc_reason     <- status$reason
  obj@meta.data$rescued       <- as.logical(status$rescue_flag)
  obj@meta.data$rescue_reason <- status$rescue_reason

  # # Enrichment files (optional)
  # ranks <- try(prep_flagged_ranks(obj, status, assay = assay), silent = TRUE)
  # if (!inherits(ranks, "try-error")) {
  #   enr <- run_enrichment(ranks, species = species, top_n_up = 200)
  #   if (!is.null(enr$ora))
  #     utils::write.csv(enr$ora,
  #                      file.path(save_dir, "removed_cells_ORA_hallmark.csv"),
  #                      row.names = FALSE)
  #   if (!is.null(enr$gsea))
  #     utils::write.csv(enr$gsea,
  #                      file.path(save_dir, "removed_cells_GSEA_hallmark.csv"),
  #                      row.names = FALSE)
  # }
  # try(removed_cells_analysis(obj,
  #                            status,
  #                            save_dir = save_dir,
  #                            species = species,
  #                            assay = assay),
  #     silent = F)

  # [6/7] Plots
  message("[6/7] Plotting QC summaries ...")
  qc_diagnostic_panel(
    obj,
    metrics,
    status_df = status,
    assay = assay,
    save_dir = save_dir,
    sample_col = sample_col,
    group_col = "auto_celltype"
  )
  qc_featuremaps(obj,
                 save_dir = save_dir,
                 prefix = "qc",
                 ncol = 4)


  if(enrichment_plots){
       validation_plots_post_annotation(
      obj,
      group.by = "auto_celltype",
      organism = species,
      assay = assay,
      slot.use = "data",
      imgDir = save_dir,
      top_n = 100,
      min_genes_for_enrich = 10,
      top_go_per_cluster = 2
    )

  }



    # [7/7] Subset & report
    message("[7/7] Subsetting object by QC & rendering report ...")

    det <- status
    stopifnot(is.data.frame(det), nrow(det) > 0)

    # Normalize to have 'qc_status'
    if (!("qc_status" %in% names(det))) {
      if ("status" %in% names(det)) {
        det$qc_status <- det$status
      } else if ("is_low_quality" %in% names(det)) {
        det$qc_status <- ifelse(det$is_low_quality, "remove", "keep")
      } else {
        stop("QC table lacks 'qc_status'/'status'/'is_low_quality'.")
      }
    }
    if (is.null(rownames(det))) stop("QC table rownames are NULL (expected cell barcodes).")

    keep_cells <- intersect(colnames(obj), rownames(det)[det$qc_status %in% c("keep","borderline")])

    if (getOption("scQCenrich.debug", FALSE)) {
      message(sprintf("[subset] keep=%d | borderline=%d | remove=%d | in_obj=%d | selected=%d",
                      sum(det$qc_status=="keep", na.rm=TRUE),
                      sum(det$qc_status=="borderline", na.rm=TRUE),
                      sum(det$qc_status=="remove", na.rm=TRUE),
                      length(intersect(colnames(obj), rownames(det))),
                      length(keep_cells)))
      if (length(keep_cells)) {
        message(sprintf("[subset] example kept: %s", paste(head(keep_cells, 5), collapse=", ")))
      }
    }

    if (length(keep_cells) == 0L) {
      utils::write.csv(det, file.path(save_dir, "qc_status.csv"), row.names = TRUE)
      stop("QC produced zero keep/borderline cells. See: ", file.path(save_dir, "qc_status.csv"))
    }

    obj_keep <- subset(obj, cells = keep_cells)

    rep_path <- if (dirname(report_file) %in% c(".", "")) file.path(save_dir, report_file) else report_file
    if (isTRUE(report_html)) {
      if (!("is_doublet" %in% colnames(obj@meta.data))) obj@meta.data$is_doublet <- FALSE
      if (!("qc_status" %in% colnames(obj@meta.data)))
        obj@meta.data$qc_status <- factor("keep", levels = c("keep","borderline","remove"))
      render_qc_report(
        obj = obj,
        metrics = metrics,
        status_df = status,
        save_dir = save_dir,
        report_file = rep_path,
        theme_bootswatch = theme_bootswatch,
        self_contained = TRUE,
        sample_col = sample_col,
        species = species
      )
    } else {
      rep_path <- NULL
    }

    invisible(list(
      obj         = obj_keep,
      obj_all     = obj,
      metrics     = metrics,
      status_df   = status,
      save_dir    = save_dir,
      report_file = if (!is.null(rep_path)) normalizePath(rep_path, winslash = "/", mustWork = FALSE) else NULL
    ))

}

# ---- helpers: doublet detection ----
detect_doublets_doubletfinder <- function(obj,
                                          assay = DefaultAssay(obj),
                                          debug = FALSE) {
  stopifnot(requireNamespace("DoubletFinder", quietly = TRUE))
  odef <- DefaultAssay(obj)
  on.exit(DefaultAssay(obj) <- odef, add = TRUE)
  DefaultAssay(obj) <- assay

  # preprocessing if needed
  if (!("pca" %in% names(obj@reductions))) {
    if (isTRUE(getOption("scQCenrich.debug")))
      message("[DF] Preprocessing (Normalize/FindVariable/Scale/PCA/Neighbors/Clusters)")
    obj <- Seurat::NormalizeData(obj, verbose = FALSE)
    obj <- Seurat::FindVariableFeatures(obj, nfeatures = 2000, verbose = FALSE)
    obj <- Seurat::ScaleData(obj, verbose = FALSE)
    obj <- Seurat::RunPCA(obj, npcs = 20, verbose = FALSE)
    obj <- Seurat::FindNeighbors(obj, dims = 1:20, verbose = FALSE)
    obj <- Seurat::FindClusters(obj, resolution = 0.5, verbose = FALSE)
  }

  sweep.list  <- DoubletFinder::paramSweep(obj, PCs = 1:20, sct = FALSE)
  sweep.stats <- DoubletFinder::summarizeSweep(sweep.list, GT = FALSE)
  bcmvn       <- DoubletFinder::find.pK(sweep.stats)
  kcol <- intersect(colnames(bcmvn), c("BCmetric", "BCmvn", "BC"))[1]
  pk   <- bcmvn$pK[which.max(bcmvn[[kcol]])]
  if (is.factor(pk))
    pk <- as.numeric(levels(pk))[pk]
  else
    pk <- as.numeric(pk)

  nExp <- max(10L, as.integer(0.075 * ncol(obj)))
  res  <- DoubletFinder::doubletFinder(
    obj,
    PCs = 1:20,
    pK = pk,
    nExp = nExp,
    pN = 0.25,
    sct = FALSE
  )

  # find the new DF columns (names contain parameters)
  md <- res@meta.data
  pan_col <- grep("^pANN_", colnames(md), value = TRUE)[1]
  cls_col <- grep("^DF\\.classifications_", colnames(md), value = TRUE)[1]

  isdbl <- setNames(grepl("(?i)doublet", md[[cls_col]]), rownames(md))
  score <- setNames(as.numeric(md[[pan_col]]), rownames(md))  # pANN in [0,1]

  out <- data.frame(
    is_doublet = isdbl,
    dbl_score  = score,
    dbl_method = "DoubletFinder",
    row.names  = names(isdbl),
    check.names = FALSE
  )
  out
}


detect_doublets_scdblfinder <- function(obj, assay = DefaultAssay(obj)) {
  stopifnot(
    requireNamespace("scDblFinder", quietly = TRUE),
    requireNamespace("SingleCellExperiment", quietly = TRUE)
  )
  counts <- Seurat::GetAssayData(obj, assay = assay, layer = "counts")
  sce <- SingleCellExperiment::SingleCellExperiment(assays = list(counts = counts))
  sce <- scDblFinder::scDblFinder(sce)
  lab   <- as.character(sce$scDblFinder.class)   # "singlet"/"doublet" etc.
  score <- as.numeric(sce$scDblFinder.score)     # probability-like
  isdbl <- setNames(grepl("(?i)doublet", lab), colnames(obj))
  out <- data.frame(
    is_doublet = isdbl,
    dbl_score  = setNames(score, colnames(obj)),
    dbl_method = "scDblFinder",
    row.names  = colnames(obj),
    check.names = FALSE
  )
  out
}
