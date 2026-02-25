# List organisms with cached MODOMICS data

Returns the names of organisms for which MODOMICS tRNA modification data
is bundled with the package. These organisms can be used with
[`modomics_mods()`](https://rnabioco.github.io/clover/reference/modomics_mods.md)
without internet access.

## Usage

``` r
modomics_organisms()
```

## Value

A character vector of organism names.

## Examples

``` r
modomics_organisms()
#> [1] "Caenorhabditis elegans"   "Drosophila melanogaster" 
#> [3] "Enterobacteria phage T4"  "Enterobacteria phage T5" 
#> [5] "Escherichia coli"         "Homo sapiens"            
#> [7] "Mus musculus"             "Saccharomyces cerevisiae"
```
