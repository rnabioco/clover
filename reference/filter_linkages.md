# Filter odds ratios for structure linkage arcs.

Convenience filter for odds ratio data that returns a tibble ready for
the `linkages` parameter of
[`plot_tRNA_structure()`](https://rnabioco.github.io/clover/reference/plot_tRNA_structure.md).
Filters by p-value, observation count, and log odds ratio magnitude,
then selects the columns needed for plotting.

## Usage

``` r
filter_linkages(data, max_p = 0.01, min_obs = 100, min_lor = 1)
```

## Arguments

- data:

  A tibble of odds ratio data, typically from
  [`clean_odds_ratios()`](https://rnabioco.github.io/clover/reference/clean_odds_ratios.md),
  with columns `pos1`, `pos2`, `log_odds_ratio`, `p_value`, and
  `total_obs`.

- max_p:

  Maximum p-value to retain. Default `0.01`.

- min_obs:

  Minimum total observations to retain. Default `100`.

- min_lor:

  Minimum absolute log odds ratio to retain. Default `1.0`.

## Value

A tibble with columns `pos1`, `pos2`, and `value` (the log odds ratio),
ready for
[`plot_tRNA_structure()`](https://rnabioco.github.io/clover/reference/plot_tRNA_structure.md).

## Examples

``` r
df <- tibble::tibble(
  pos1 = c(20, 34, 10),
  pos2 = c(34, 58, 45),
  odds_ratio = c(4.0, 0.3, 1.1),
  log_odds_ratio = c(1.4, -1.2, 0.1),
  p_value = c(0.001, 0.005, 0.5),
  total_obs = c(200, 150, 50)
)
filter_linkages(df)
#> # A tibble: 2 × 3
#>    pos1  pos2 value
#>   <dbl> <dbl> <dbl>
#> 1    20    34   1.4
#> 2    34    58  -1.2
```
