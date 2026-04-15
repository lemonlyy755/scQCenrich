# List the tissue names supported by the PanglaoDB marker database

Prints (and invisibly returns) the character vector of tissue names that
can be passed to the `tissue` argument of
[`run_qc_pipeline`](https://lemonlyy755.github.io/scQCenrich/reference/run_qc_pipeline.md)
or `panglao_signatures`. Supplying the correct tissue name(s) narrows
the signature lookup to cell types relevant to your sample and
substantially improves annotation accuracy.

## Usage

``` r
list_panglao_tissues()
```

## Value

A character vector of supported tissue names (invisibly).

## Examples

``` r
list_panglao_tissues()
#> Supported PanglaoDB tissue names (pass to `tissue` argument):
#>   1. Adrenal glands
#>   2. Blood
#>   3. Bone
#>   4. Brain
#>   5. Connective tissue
#>   6. Embryo
#>   7. Epithelium
#>   8. Eye
#>   9. GI tract
#>   10. Heart
#>   11. Immune system
#>   12. Kidney
#>   13. Liver
#>   14. Lungs
#>   15. Mammary gland
#>   16. Olfactory system
#>   17. Oral cavity
#>   18. Pancreas
#>   19. Parathyroid glands
#>   20. Placenta
#>   21. Reproductive
#>   22. Skeletal muscle
#>   23. Skin
#>   24. Smooth muscle
#>   25. Thymus
#>   26. Thyroid
#>   27. Urinary bladder
#>   28. Vasculature
#>   29. Zygote
#> 
#> Example: run_qc_pipeline(..., tissue = c("Blood", "Immune system"))
```
