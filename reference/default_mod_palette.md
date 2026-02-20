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
#>       m1A       m1G       m2A       m2G     m2,2G       m3C       m5C       m5U 
#> "#E41A1C" "#984EA3" "#E68A00" "#66C2A5" "#8DA0CB" "#B3DE69" "#377EB8" "#FF7F00" 
#>       m6A     m6t6A       m7G         D         Y        Ym         I       t6A 
#> "#FB9A99" "#D95F02" "#4DAF4A" "#FFFF33" "#A65628" "#C49A6C" "#F781BF" "#FC8D62" 
#>         Q        Cm        Gm        Um        Am      ac4C     acp3U       s2C 
#> "#1B9E77" "#E78AC3" "#A6D854" "#FFD92F" "#E5C494" "#B2DF8A" "#CAB2D6" "#33A02C" 
#>       s4U       i6A    ms2i6A     mnm5U   mnm5s2U     cmo5U  cmnm5s2U   cmnm5Um 
#> "#6A3D9A" "#FF6699" "#B15928" "#7570B3" "#E7298A" "#66A61E" "#A6761D" "#D4A76A" 
```
