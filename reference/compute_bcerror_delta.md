# Compute per-position delta between two conditions.

Take a summarized bcerror tibble (with columns for reference, position,
condition, and a value such as mean error rate) and compute the
difference between two conditions at each position.

## Usage

``` r
compute_bcerror_delta(
  data,
  delta,
  value_col = "mean_error",
  condition_col = "condition"
)
```

## Arguments

- data:

  A summarized bcerror tibble with at least `ref`, `pos`, and columns
  named by `value_col` and `condition_col`.

- delta:

  A bare expression of the form `lhs - rhs`, where `lhs` and `rhs` are
  condition levels. The result is `lhs - rhs`.

- value_col:

  Column name (string) containing the values to pivot. Default
  `"mean_error"`.

- condition_col:

  Column name (string) containing condition labels. Default
  `"condition"`.

## Value

A tibble with columns `ref`, `pos`, one column per condition level, and
`delta` (the computed difference).

## Examples

``` r
df <- tidyr::expand_grid(
  ref = c("tRNA-Ala", "tRNA-Gly"),
  pos = 1:5,
  condition = c("wt", "mut")
)
df$mean_error <- runif(nrow(df), 0, 0.3)
compute_bcerror_delta(df, delta = wt - mut)
#> # A tibble: 10 × 5
#>    ref        pos      wt    mut    delta
#>    <chr>    <int>   <dbl>  <dbl>    <dbl>
#>  1 tRNA-Ala     1 0.158   0.180  -0.0218 
#>  2 tRNA-Ala     2 0.0784  0.0870 -0.00860
#>  3 tRNA-Ala     3 0.144   0.276  -0.132  
#>  4 tRNA-Ala     4 0.120   0.0640  0.0563 
#>  5 tRNA-Ala     5 0.202   0.0176  0.184  
#>  6 tRNA-Gly     1 0.299   0.0447  0.254  
#>  7 tRNA-Gly     2 0.156   0.254  -0.0983 
#>  8 tRNA-Gly     3 0.215   0.0724  0.143  
#>  9 tRNA-Gly     4 0.164   0.250  -0.0863 
#> 10 tRNA-Gly     5 0.00839 0.141  -0.132  
```
