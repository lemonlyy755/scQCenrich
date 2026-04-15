# Benchmark Results

This file summarizes benchmark artifacts currently present under `benchmarking`. It does not rerun benchmarks, infer results from exploratory outputs, or generate new benchmark data.

| Dataset | Input shape | Benchmark artifact(s) found | Status | Key takeaway |
| --- | --- | --- | --- | --- |
| brain | `filtered_feature_bc_matrix/` + `testloom.rds` | None | Missing summary artifact | No benchmark summary artifact found in current repo. |
| GSE131907_Lung_Cancer | Annotation table, raw UMI matrix, normalized log2TPM matrix, feature summary spreadsheet | `GSE131907_2sample_benchmark_results.rds`, `GSE131907_4sample_benchmark_results.rds`, `GSE131907_Benchmark_Report.pdf`, `GSE131907_4S_Benchmark_Report.pdf`, `GSE131907_Retention_BarPlot.png`, `GSE131907_Concordance_UpSet.png`, `GSE131907_4S_Retention_BarPlot.png`, `GSE131907_4S_CellType_Retention.png`, `GSE131907_4S_Concordance_UpSet.png` | Summarized | Both benchmark runs are available; `miQC` or `SampleQC` collapses to zero in at least one run. |
| GSE46980_CombinedMoleculeCounts.tab | Plate-based counts table | `GSE46980_benchmark_results.txt`, `GSE46980_benchmark_eval_table.csv`, `QC_Benchmarking_GSE46980.png` | Summarized | `scater_MAD3` has the highest recall, while `Standard_MT5_nF200` has the best precision and accuracy in the recorded summary. |
| GSE75790 | Empty directory | None | Missing summary artifact | No benchmark summary artifact found in current repo. |
| GSE92495 | Empty directory | None | Missing summary artifact | No benchmark summary artifact found in current repo. |
| heart | `filtered_feature_bc_matrix/` + `testloom.rds` | None | Missing summary artifact | No benchmark summary artifact found in current repo. |
| lungCancer7k | `filtered_feature_bc_matrix/` + `testloom.rds` | `ground_truth_benchmark_results.txt` | Summarized | The recorded ground-truth summary reports that scQCenrich retained more annotated intact cells than the standard threshold. |
| melanoma10k | `filtered_feature_bc_matrix/` + `testloom.rds` | None | Missing summary artifact | No benchmark summary artifact found in current repo. |
| overian | `filtered_feature_bc_matrix/` + `testloom.rds` | None | Missing summary artifact | No benchmark summary artifact found in current repo. |
| pbmc | Filtered 10x matrix, raw 10x matrix, `pbmc1k_v3_cr.loom`, `testloom.rds` | None | Missing summary artifact | No benchmark summary artifact found in current repo. |
| pbmc10k | `filtered_feature_bc_matrix/` + `testloom.rds` | None | Missing summary artifact | No benchmark summary artifact found in current repo. |
| pbmc3k | Legacy `filtered_gene_bc_matrices/` directory | None | Missing summary artifact | No benchmark summary artifact found in current repo. |

## brain

- Input shape: `filtered_feature_bc_matrix/` and `testloom.rds`.
- Benchmark artifact(s) found: none.
- Status: No benchmark summary artifact found in current repo.

## GSE131907_Lung_Cancer

- Input shape: `GSE131907_Lung_Cancer_cell_annotation.txt`, `GSE131907_Lung_Cancer_Feature_Summary.xlsx`, `GSE131907_Lung_Cancer_raw_UMI_matrix.txt`, `GSE131907_Lung_Cancer_raw_UMI_matrix.rds`, `GSE131907_Lung_Cancer_normalized_log2TPM_matrix.txt`, and `GSE131907_Lung_Cancer_normalized_log2TPM_matrix.rds`.
- Benchmark artifact(s) found: `GSE131907_2sample_benchmark_results.rds`, `GSE131907_4sample_benchmark_results.rds`, `GSE131907_Benchmark_Report.pdf`, `GSE131907_4S_Benchmark_Report.pdf`, `GSE131907_Retention_BarPlot.png`, `GSE131907_Concordance_UpSet.png`, `GSE131907_4S_Retention_BarPlot.png`, `GSE131907_4S_CellType_Retention.png`, `GSE131907_4S_Concordance_UpSet.png`.

### 2-sample run

Per-sample retention percentages from `GSE131907_2sample_benchmark_results.rds`:

| Method | Sample | Retention (%) |
| --- | --- | ---: |
| scQCenrich | LUNG_N20 | 87.96 |
| scQCenrich | LUNG_T31 | 92.11 |
| SampleQC | LUNG_N20 | 0.00 |
| SampleQC | LUNG_T31 | 0.00 |
| Standard_MT5 | LUNG_N20 | 45.24 |
| Standard_MT5 | LUNG_T31 | 86.95 |

- Total rows in `pass_logic`: 11,376.
- Overall kept counts from `colSums(pass_logic)`: `scQCenrich = 10238`, `SampleQC = 0`, `Standard_MT5 = 7473`.
- Recorded takeaway: scQCenrich keeps most cells in both samples in this stored two-sample benchmark, while SampleQC keeps none.

### 4-sample run

Per-sample totals from `GSE131907_4sample_benchmark_results.rds`:

| Method | Sample | Total | Kept | Retention (%) |
| --- | --- | ---: | ---: | ---: |
| scQCenrich | LUNG_N18 | 4628 | 4298 | 92.87 |
| scQCenrich | LUNG_N20 | 5798 | 5250 | 90.55 |
| scQCenrich | LUNG_T18 | 3705 | 2962 | 79.95 |
| scQCenrich | LUNG_T31 | 5578 | 5228 | 93.73 |
| SampleQC | LUNG_N18 | 4628 | 4009 | 86.62 |
| SampleQC | LUNG_N20 | 5798 | 5457 | 94.12 |
| SampleQC | LUNG_T18 | 3705 | 2383 | 64.32 |
| SampleQC | LUNG_T31 | 5578 | 4271 | 76.57 |
| scater | LUNG_N18 | 4628 | 4407 | 95.22 |
| scater | LUNG_N20 | 5798 | 5650 | 97.45 |
| scater | LUNG_T18 | 3705 | 3116 | 84.10 |
| scater | LUNG_T31 | 5578 | 5176 | 92.79 |
| miQC | LUNG_N18 | 4628 | 0 | 0.00 |
| miQC | LUNG_N20 | 5798 | 0 | 0.00 |
| miQC | LUNG_T18 | 3705 | 0 | 0.00 |
| miQC | LUNG_T31 | 5578 | 0 | 0.00 |
| Standard_MT10 | LUNG_N18 | 4628 | 4448 | 96.11 |
| Standard_MT10 | LUNG_N20 | 5798 | 5513 | 95.08 |
| Standard_MT10 | LUNG_T18 | 3705 | 3025 | 81.65 |
| Standard_MT10 | LUNG_T31 | 5578 | 5455 | 97.79 |

- Total rows in `pass_logic`: 19,709.
- Overall kept counts from `colSums(pass_logic)`: `scQCenrich = 17738`, `SampleQC = 16120`, `scater = 18349`, `miQC = 0`, `Standard_MT10 = 18441`.
- Recorded takeaway: the four-sample artifact shows broad retention across methods except `miQC`, which keeps zero cells in this stored run.

## GSE46980_CombinedMoleculeCounts.tab

- Input shape: `GSE46980_CombinedMoleculeCounts.tab`.
- Benchmark artifact(s) found: `GSE46980_benchmark_results.txt`, `GSE46980_benchmark_eval_table.csv`, `QC_Benchmarking_GSE46980.png`.
- Ground-truth totals from `GSE46980_benchmark_results.txt`: 50 annotated intact cells and 46 annotated dead/empty cells.

| Method | Kept Total | TP | FP | FN | Recall (%) | Precision (%) | F1-score | Accuracy (%) |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Standard_MT5_nF200 | 23 | 23 | 0 | 27 | 46.00 | 100.00 | 0.630 | 71.88 |
| scater_MAD3 | 96 | 50 | 46 | 0 | 100.00 | 52.08 | 0.685 | 52.08 |
| scQCenrich | 86 | 43 | 43 | 7 | 86.00 | 50.00 | 0.632 | 47.92 |

- `GSE46980_benchmark_results.txt` also records that `miQC` failed and kept 0 cells.
- `GSE46980_benchmark_results.txt` also records that `SampleQC` failed and kept 0 cells.
- Recorded takeaway: `scater_MAD3` maximizes recall, while `Standard_MT5_nF200` is the most conservative method and the only one with zero false positives in the saved evaluation table.

## GSE75790

- Input shape: empty directory.
- Benchmark artifact(s) found: none.
- Status: No benchmark summary artifact found in current repo.

## GSE92495

- Input shape: empty directory.
- Benchmark artifact(s) found: none.
- Status: No benchmark summary artifact found in current repo.

## heart

- Input shape: `filtered_feature_bc_matrix/` and `testloom.rds`.
- Benchmark artifact(s) found: none.
- Status: No benchmark summary artifact found in current repo.

## lungCancer7k

- Input shape: `filtered_feature_bc_matrix/` and `testloom.rds`.
- Benchmark artifact(s) found: `ground_truth_benchmark_results.txt`.

Recorded totals from `ground_truth_benchmark_results.txt`:

| Metric | Value |
| --- | --- |
| Total annotated intact cells | 7092 |
| Standard threshold kept | 5162 (72.79%) |
| scQCenrich kept | 6406 (90.33%) |
| Rescued intact cells vs standard threshold | 1244 |

- Recorded takeaway: the saved ground-truth summary attributes 1,244 additional intact-cell rescues to scQCenrich relative to the standard mitochondrial threshold.

## melanoma10k

- Input shape: `filtered_feature_bc_matrix/` and `testloom.rds`.
- Benchmark artifact(s) found: none.
- Status: No benchmark summary artifact found in current repo.

## overian

- Input shape: `filtered_feature_bc_matrix/` and `testloom.rds`.
- Benchmark artifact(s) found: none.
- Status: No benchmark summary artifact found in current repo.

## pbmc

- Input shape: `filtered_feature_bc_matrix/`, `raw_feature_bc_matrix/`, `pbmc1k_v3_cr.loom`, and `testloom.rds`.
- Benchmark artifact(s) found: none.
- Status: No benchmark summary artifact found in current repo.

## pbmc10k

- Input shape: `filtered_feature_bc_matrix/` and `testloom.rds`.
- Benchmark artifact(s) found: none.
- Status: No benchmark summary artifact found in current repo.

## pbmc3k

- Input shape: `filtered_gene_bc_matrices/`.
- Benchmark artifact(s) found: none.
- Status: No benchmark summary artifact found in current repo.

## Appendix: Standalone Benchmark-Only Artifacts

### E-MTAB-2805

- Artifact: `emtab2805_fluidigm_benchmark.txt`.
- Benchmark framing: false-positive analysis on verified intact cells from the Buettner et al. mESC cell-cycle dataset.
- Total verified cells recorded: 297.

| Method | Kept | Lost intact cells |
| --- | ---: | ---: |
| Hard Threshold (MT < 5%) | 297 (100.00%) | 0 |
| Hard Threshold (MT < 10%) | 297 (100.00%) | 0 |
| MAD Outlier Filtering | 297 (100.00%) | 0 |

- Recorded takeaway: the saved benchmark reports no false-positive removal for any of the listed methods on this fully verified intact-cell dataset.
