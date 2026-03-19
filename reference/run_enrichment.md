# Hallmark enrichment (ORA + GSEA) for dropped cells

Hallmark enrichment (ORA + GSEA) for dropped cells

## Usage

``` r
run_enrichment(ranks, species = c("mouse", "human"), top_n_up = 200)
```

## Arguments

- ranks:

  Named numeric vector of gene statistics

- species:

  'mouse' or 'human'

- top_n_up:

  Number of top up-genes to use for ORA

## Value

List with data.frames: \$ora and \$gsea (or NULLs if unavailable)
