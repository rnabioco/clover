# List available tRNA structures for an organism

Returns the names of tRNAs for which cloverleaf structure SVGs are
bundled with the package for the given organism.

## Usage

``` r
structure_trnas(organism)
```

## Arguments

- organism:

  Character string specifying the organism name (e.g.,
  `"Escherichia coli"`). Use
  [`structure_organisms()`](https://rnabioco.github.io/clover/reference/structure_organisms.md)
  to list available organisms.

## Value

A character vector of tRNA names.

## Examples

``` r
structure_trnas("Escherichia coli")
#>  [1] "tRNA-Ala-GGC"  "tRNA-Ala-TGC"  "tRNA-Arg-ACG"  "tRNA-Arg-CCG" 
#>  [5] "tRNA-Arg-CCT"  "tRNA-Arg-TCT"  "tRNA-Asn-GTT"  "tRNA-Asp-GTC" 
#>  [9] "tRNA-Cys-GCA"  "tRNA-Gln-CTG"  "tRNA-Gln-TTG"  "tRNA-Glu-TTC" 
#> [13] "tRNA-Gly-CCC"  "tRNA-Gly-GCC"  "tRNA-Gly-TCC"  "tRNA-His-GTG" 
#> [17] "tRNA-Ile-GAT"  "tRNA-Ile2-CAT" "tRNA-Leu-CAA"  "tRNA-Leu-CAG" 
#> [21] "tRNA-Leu-GAG"  "tRNA-Leu-TAA"  "tRNA-Leu-TAG"  "tRNA-Lys-TTT" 
#> [25] "tRNA-Met-CAT"  "tRNA-Phe-GAA"  "tRNA-Pro-CGG"  "tRNA-Pro-GGG" 
#> [29] "tRNA-Pro-TGG"  "tRNA-SeC-TCA"  "tRNA-Ser-CGA"  "tRNA-Ser-GCT" 
#> [33] "tRNA-Ser-GGA"  "tRNA-Ser-TGA"  "tRNA-Thr-CGT"  "tRNA-Thr-GGT" 
#> [37] "tRNA-Thr-TGT"  "tRNA-Trp-CCA"  "tRNA-Tyr-GTA"  "tRNA-Val-GAC" 
#> [41] "tRNA-Val-TAC"  "tRNA-fMet-CAT"
```
