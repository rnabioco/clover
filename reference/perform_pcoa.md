# Perform PCoA on a rewiring matrix.

Run classical multidimensional scaling (PCoA) on a Euclidean distance
matrix derived from the ROR matrix.

## Usage

``` r
perform_pcoa(mat, k = 2)
```

## Arguments

- mat:

  A numeric matrix from
  [`prepare_rewiring_matrix()`](https://rnabioco.github.io/clover/reference/prepare_rewiring_matrix.md).

- k:

  Number of dimensions. Default `2`.

## Value

A list with elements:

- `coordinates`: a tibble with `isodecoder`, `PC1`, `PC2`, ...

- `variance_explained`: numeric vector of percent variance explained for
  each dimension

- `eigenvalues`: full eigenvalue vector from
  [`stats::cmdscale()`](https://rdrr.io/r/stats/cmdscale.html)

## Examples

``` r
mat <- matrix(
  c(1.5, -0.8, 0.3, 2.1, 0.5, -1.2),
  nrow = 3,
  dimnames = list(
    c("tRNA-Ala", "tRNA-Gly", "tRNA-Ser"),
    c("20_vs_34", "34_vs_58")
  )
)
perform_pcoa(mat)
#> $coordinates
#> # A tibble: 3 × 3
#>      PC1    PC2 isodecoder
#>    <dbl>  <dbl> <chr>     
#> 1  1.98   0.299 tRNA-Ala  
#> 2 -0.484 -1.03  tRNA-Gly  
#> 3 -1.50   0.726 tRNA-Ser  
#> 
#> $variance_explained
#> [1] 79.39345 20.60655
#> 
#> $eigenvalues
#> [1]  6.425576e+00  1.667757e+00 -2.664535e-15
#> 
```
