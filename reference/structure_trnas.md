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
#>  [1] "tRNA-Ala-GGC" "tRNA-Ala-TGC" "tRNA-Arg-ACG" "tRNA-Arg-CCG" "tRNA-Arg-CCT"
#>  [6] "tRNA-Arg-TCT" "tRNA-Asn-GTT" "tRNA-Asp-GTC" "tRNA-Cys-GCA" "tRNA-Gln-CTG"
#> [11] "tRNA-Gln-TTG" "tRNA-Glu-TTC" "tRNA-Gly-CCC" "tRNA-Gly-GCC" "tRNA-Gly-TCC"
#> [16] "tRNA-His-GTG" "tRNA-Ile-GAT" "tRNA-Leu-CAA" "tRNA-Leu-CAG" "tRNA-Leu-GAG"
#> [21] "tRNA-Leu-TAA" "tRNA-Leu-TAG" "tRNA-Lys-TTT" "tRNA-Met-CAT" "tRNA-Phe-GAA"
#> [26] "tRNA-Pro-CGG" "tRNA-Pro-GGG" "tRNA-Pro-TGG" "tRNA-Ser-CGA" "tRNA-Ser-GCT"
#> [31] "tRNA-Ser-GGA" "tRNA-Ser-TGA" "tRNA-Thr-CGT" "tRNA-Thr-GGT" "tRNA-Thr-TGT"
#> [36] "tRNA-Trp-CCA" "tRNA-Tyr-GTA" "tRNA-Val-GAC" "tRNA-Val-TAC"
```
