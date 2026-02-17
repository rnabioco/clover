# Prepare a matrix for rewiring PCoA analysis.

Build a wide matrix from isodecoder-level relative odds ratios suitable
for PCoA or distance calculations. Rows are isodecoders, columns are
position pair comparisons.

## Usage

``` r
prepare_rewiring_matrix(data, sig_only = TRUE, value_cap = 10)
```

## Arguments

- data:

  A tibble with columns `isodecoder`, `pos1`, `pos2`, and `ror`.
  Typically the output of
  [`compute_ror_isodecoder()`](https://rnabioco.github.io/clover/reference/compute_ror_isodecoder.md).

- sig_only:

  Logical; filter to rows where `significant` is `TRUE`? Default `TRUE`.

- value_cap:

  Maximum absolute ROR value. Values beyond this are capped. Default
  `10`.

## Value

A numeric matrix with isodecoders as row names and position pair labels
(`pos1_vs_pos2`) as column names.

## Examples

``` r
df <- tibble::tibble(
  isodecoder = rep(c("tRNA-Ala", "tRNA-Gly"), each = 2),
  pos1 = c(20, 34, 20, 34),
  pos2 = c(34, 58, 34, 58),
  ror = c(1.5, -0.8, 0.3, 2.1),
  significant = c(TRUE, TRUE, TRUE, TRUE)
)
prepare_rewiring_matrix(df)
#>          20_vs_34 34_vs_58
#> tRNA-Ala      1.5     -0.8
#> tRNA-Gly      0.3      2.1
```
