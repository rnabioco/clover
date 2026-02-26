# List organisms with bundled tRNA structure SVGs

Returns the names of organisms for which tRNA cloverleaf structure SVGs
are bundled with the package. These organisms can be used with
[`plot_tRNA_structure()`](https://rnabioco.github.io/clover/reference/plot_tRNA_structure.md).

## Usage

``` r
structure_organisms()
```

## Value

A character vector of organism names.

## Examples

``` r
structure_organisms()
#> [1] "Escherichia coli"         "Homo sapiens"            
#> [3] "Saccharomyces cerevisiae" "T4 phage"                
```
