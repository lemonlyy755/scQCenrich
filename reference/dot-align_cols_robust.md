# Try to align matrix columns to target cells with multiple guarded passes Returns a matrix subset/reordered to target_cells; attributes include alignment stats and the pass used. Never introduces regressions: if a pass produces duplicate/ambiguous matches, it is skipped.

Try to align matrix columns to target cells with multiple guarded passes
Returns a matrix subset/reordered to target_cells; attributes include
alignment stats and the pass used. Never introduces regressions: if a
pass produces duplicate/ambiguous matches, it is skipped.

## Usage

``` r
.align_cols_robust(mat, target_cells, tag = "ext")
```
