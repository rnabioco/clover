# Clean odds ratio data.

Prepare odds ratio data for downstream analysis by capping infinite log
odds ratio values and optionally removing adapter positions.

## Usage

``` r
clean_odds_ratios(data, cap_inf = 2)
```

## Arguments

- data:

  A tibble of odds ratio data with columns `ref`, `odds_ratio`, and
  `log_odds_ratio`.

- cap_inf:

  Numeric value added to (or subtracted from) the finite range when
  capping positive (or negative) infinite log odds ratios. Default `2`.

## Value

A tibble with a new `log_or_clean` column where infinite values have
been capped.

## Examples

``` r
df <- tibble::tibble(
  ref = "tRNA-Ala-AGC-1-1",
  pos1 = c(20, 34),
  pos2 = c(34, 58),
  odds_ratio = c(2.5, Inf),
  log_odds_ratio = c(0.92, Inf),
  p_value = c(0.001, 0.01),
  total_obs = c(200, 150)
)
clean_odds_ratios(df)
#> # A tibble: 2 × 8
#>   ref        pos1  pos2 odds_ratio log_odds_ratio p_value total_obs log_or_clean
#>   <chr>     <dbl> <dbl>      <dbl>          <dbl>   <dbl>     <dbl>        <dbl>
#> 1 tRNA-Ala…    20    34        2.5           0.92   0.001       200         0.92
#> 2 tRNA-Ala…    34    58      Inf           Inf      0.01        150         2.92
```
