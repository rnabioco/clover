# Calculate rewiring scores from a ROR matrix.

Summarize the magnitude and extent of rewiring for each isodecoder using
Euclidean magnitude, mean and max absolute change, and the number of
non-zero comparisons.

## Usage

``` r
calculate_rewiring_scores(mat)
```

## Arguments

- mat:

  A numeric matrix from
  [`prepare_rewiring_matrix()`](https://rnabioco.github.io/clover/reference/prepare_rewiring_matrix.md),
  with isodecoders as row names.

## Value

A tibble with columns: `isodecoder`, `euclidean_magnitude`,
`mean_abs_change`, `max_abs_change`, and `n_nonzero`, sorted by
descending `euclidean_magnitude`.

## Examples

``` r
mat <- matrix(
  c(1.5, -0.8, 0.3, 2.1),
  nrow = 2,
  dimnames = list(c("tRNA-Ala", "tRNA-Gly"), c("20_vs_34", "34_vs_58"))
)
calculate_rewiring_scores(mat)
#> # A tibble: 2 × 5
#>   isodecoder euclidean_magnitude mean_abs_change max_abs_change n_nonzero
#>   <chr>                    <dbl>           <dbl>          <dbl>     <dbl>
#> 1 tRNA-Gly                  2.25            1.45            2.1         2
#> 2 tRNA-Ala                  1.53            0.9             1.5         2
```
