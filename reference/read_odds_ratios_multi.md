# Read odds ratio files for multiple samples.

Read and combine odds ratio files from multiple samples into a single
tibble with a `sample_id` column.

## Usage

``` r
read_odds_ratios_multi(paths)
```

## Arguments

- paths:

  A named character vector of file paths. Names are used as sample
  identifiers.

## Value

A tibble with all samples combined and a `sample_id` column.

## Examples

``` r
paths <- c(
  "wt-15-ctl-01" = clover_example(
    "ecoli/summary/tables/wt-15-ctl-01/wt-15-ctl-01.odds_ratios.tsv.gz"
  ),
  "wt-15-ctl-02" = clover_example(
    "ecoli/summary/tables/wt-15-ctl-02/wt-15-ctl-02.odds_ratios.tsv.gz"
  )
)
read_odds_ratios_multi(paths)
#> # A tibble: 22,192 × 8
#>    sample_id    ref       pos1  pos2 odds_ratio log_odds_ratio p_value total_obs
#>    <chr>        <chr>    <dbl> <dbl>      <dbl>          <dbl>   <dbl>     <dbl>
#>  1 wt-15-ctl-01 host-tR…    20    23          0          -23.0       1      1488
#>  2 wt-15-ctl-01 host-tR…    20    26          0          -23.0       1      1488
#>  3 wt-15-ctl-01 host-tR…    20    28          0          -23.0       1      1488
#>  4 wt-15-ctl-01 host-tR…    20    31          0          -23.0       1      1488
#>  5 wt-15-ctl-01 host-tR…    20    32          0          -23.0       1      1488
#>  6 wt-15-ctl-01 host-tR…    20    35          0          -23.0       1      1488
#>  7 wt-15-ctl-01 host-tR…    20    36          0          -23.0       1      1488
#>  8 wt-15-ctl-01 host-tR…    20    37          0          -23.0       1      1488
#>  9 wt-15-ctl-01 host-tR…    20    39          0          -23.0       1      1488
#> 10 wt-15-ctl-01 host-tR…    20    40          0          -23.0       1      1488
#> # ℹ 22,182 more rows
```
