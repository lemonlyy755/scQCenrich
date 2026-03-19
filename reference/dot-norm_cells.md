# Normalize cell barcodes in several safe modes (used for cross-obj alignment)

Normalize cell barcodes in several safe modes (used for cross-obj
alignment)

## Usage

``` r
.norm_cells(cells, mode = c("basic", "core16"))
```

## Arguments

- cells:

  character vector of cell names

- mode:

  "basic" strips sample prefix (before ':'), strips trailing 'x', and
  strips trailing '-' (e.g., '-1'); "core16" returns just the last 16
  A/C/G/T characters if present.
