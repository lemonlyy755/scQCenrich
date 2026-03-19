# Package index

## All functions

- [`applyQCfilter()`](https://lemonlyy755.github.io/scQCenrich/reference/applyQCfilter.md)
  : Subset a Seurat object by QC status (Seurat v5-safe)
- [`apply_tuned_clustering()`](https://lemonlyy755.github.io/scQCenrich/reference/apply_tuned_clustering.md)
  : Apply chosen k/res to an object (returns updated object)
- [`auto_annotate()`](https://lemonlyy755.github.io/scQCenrich/reference/auto_annotate.md)
  : Auto-annotate cell types (multi-backend wrapper)
- [`calcQCmetrics()`](https://lemonlyy755.github.io/scQCenrich/reference/calcQCmetrics.md)
  : Calculate QC metrics (UMIs, genes, mito%, MALAT1, stress; optional
  splicing metrics) Computes nCount/nFeature, pctMT, MALAT1 fraction,
  stress_score (case-insensitive SYMBOL match), and optional splicing
  metrics from same object layers or external Seurat objects. Compatible
  with Seurat v5 layers and v4 slots. Returns a data.frame and (if
  add_to_meta=TRUE) writes \*\_QC fields into obj@meta.data (no
  regression in field names).
- [`default_stress_genes()`](https://lemonlyy755.github.io/scQCenrich/reference/default_stress_genes.md)
  : Default stress genes
- [`flagLowQuality()`](https://lemonlyy755.github.io/scQCenrich/reference/flagLowQuality.md)
  : Flag low-quality cells
- [`generic_signatures()`](https://lemonlyy755.github.io/scQCenrich/reference/generic_signatures.md)
  : Generic, cross-tissue marker signatures (broad/extended)
- [`prep_flagged_ranks()`](https://lemonlyy755.github.io/scQCenrich/reference/prep_flagged_ranks.md)
  : Prepare ranked gene stats for dropped vs kept (logCPM mean diff)
- [`qc_featuremaps()`](https://lemonlyy755.github.io/scQCenrich/reference/qc_featuremaps.md)
  : Make & save UMAP featuremaps for QC metrics
- [`removed_cells_analysis()`](https://lemonlyy755.github.io/scQCenrich/reference/removed_cells_analysis.md)
  : Removed cells clustering + enrichment (used by report)
- [`render_qc_report()`](https://lemonlyy755.github.io/scQCenrich/reference/render_qc_report.md)
  : Render a pretty HTML QC report
- [`rescue_by_coherence()`](https://lemonlyy755.github.io/scQCenrich/reference/rescue_by_coherence.md)
  : Rescue borderline/removed cells by cluster coherence
- [`run_enrichment()`](https://lemonlyy755.github.io/scQCenrich/reference/run_enrichment.md)
  : Hallmark enrichment (ORA + GSEA) for dropped cells
- [`run_qc_pipeline()`](https://lemonlyy755.github.io/scQCenrich/reference/run_qc_pipeline.md)
  : High-level wrapper that computes QC metrics, flags low-quality
  cells, (optionally) detects doublets and performs lightweight
  annotation, and renders an HTML report.
- [`scq_plot_auto_confidence()`](https://lemonlyy755.github.io/scQCenrich/reference/scq_plot_auto_confidence.md)
  : FeaturePlot-like map of auto-annotation confidence
- [`toy_seu`](https://lemonlyy755.github.io/scQCenrich/reference/toy_seu.md)
  : toy_seu
- [`tune_clustering()`](https://lemonlyy755.github.io/scQCenrich/reference/tune_clustering.md)
  : Grid-search to preserve small clusters (k.param & resolution)
