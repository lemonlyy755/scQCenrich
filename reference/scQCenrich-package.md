# scQCenrich: State-aware QC for single-cell RNA-seq (with spliced/unspliced support)

Features:

- QC outlier detection (GMM/HDBSCAN/threshold) over mito%, genes,
  stress, MALAT1, and (optionally) spliced/unspliced metrics
  (intronic_frac, unspliced/spliced fractions & ratio).

- Rescue of coherent, non-extreme flagged cells (3-tier
  keep/borderline/drop).

- Doublet detection (DoubletFinder / scDblFinder).

- Automatic cell-type annotation

- Presentable HTML report (QC×label, removed-cell UMAP, enrichment).

## Details

See README for how to set `spliced` / `unspliced` inputs (assay or
assay:layer).

## Author

**Maintainer**: Yuanyuan Liu <18211007107@163.com>
