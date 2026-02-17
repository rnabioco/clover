# Aggregate odds ratios to isodecoder level.

Collapse per-gene odds ratios to isodecoder level by removing gene copy
numbers from the reference name and summarizing across copies.

## Usage

``` r
aggregate_or_isodecoder(data, pattern = "-\\d+-\\d+$")
```

## Arguments

- data:

  A tibble with columns `ref`, `pos1`, `pos2`, `log_or_clean` (from
  [`clean_odds_ratios()`](https://rnabioco.github.io/clover/reference/clean_odds_ratios.md)),
  `p_value`, and `total_obs`.

- pattern:

  Regular expression used to strip gene copy numbers from `ref`. Default
  `"-\\d+-\\d+$"` matches trailing `-<copy>-<gene>` suffixes.

## Value

A tibble with columns: `isodecoder`, `pos1`, `pos2`, `mean_or`,
`mean_log_or`, `sd_log_or`, `min_pval`, `total_reads`, and `n_copies`.

## Examples

``` r
df <- tibble::tibble(
  ref = c("tRNA-Ala-AGC-1-1", "tRNA-Ala-AGC-2-1"),
  pos1 = c(20, 20),
  pos2 = c(34, 34),
  odds_ratio = c(2.5, 3.0),
  log_or_clean = c(0.92, 1.10),
  p_value = c(0.001, 0.01),
  total_obs = c(200, 150)
)
aggregate_or_isodecoder(df)
#> # A tibble: 1 × 9
#>   isodecoder    pos1  pos2 mean_or mean_log_or sd_log_or min_pval total_reads
#>   <chr>        <dbl> <dbl>   <dbl>       <dbl>     <dbl>    <dbl>       <dbl>
#> 1 tRNA-Ala-AGC    20    34    2.75        1.01     0.127    0.001         350
#> # ℹ 1 more variable: n_copies <int>
```
