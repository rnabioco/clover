# Default modification color palette

Returns a named character vector of colors for common tRNA
modifications, suitable for use with
[`plot_tRNA_structure()`](https://rnabioco.github.io/clover/reference/plot_tRNA_structure.md).

## Usage

``` r
default_mod_palette()
```

## Value

A named character vector mapping modification names to hex colors.

## Examples

``` r
default_mod_palette()
#>       m1A       m5C       m7G       m1G       m5U         D   pseudoU         I 
#> "#E41A1C" "#377EB8" "#4DAF4A" "#984EA3" "#FF7F00" "#FFFF33" "#A65628" "#F781BF" 
#>       m2G       t6A      m22G        Cm        Gm        Um        Am 
#> "#66C2A5" "#FC8D62" "#8DA0CB" "#E78AC3" "#A6D854" "#FFD92F" "#E5C494" 
```
