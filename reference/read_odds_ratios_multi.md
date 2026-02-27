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
or_file1 <- paste0(
  "ecoli/summary/tables/wt-15-ctl-01/",
  "wt-15-ctl-01.odds_ratios_filtered.tsv.gz"
)
or_file2 <- paste0(
  "ecoli/summary/tables/wt-15-ctl-02/",
  "wt-15-ctl-02.odds_ratios_filtered.tsv.gz"
)
paths <- c(
  "wt-15-ctl-01" = clover_example(or_file1),
  "wt-15-ctl-02" = clover_example(or_file2)
)
read_odds_ratios_multi(paths)
#> # A tibble: 1,585 × 17
#>    sample_id    ref     pos1  pos2    n00   n01   n10   n11 total_obs odds_ratio
#>    <chr>        <chr>   <chr> <chr> <dbl> <dbl> <dbl> <dbl>     <dbl>      <dbl>
#>  1 wt-15-ctl-01 host-t… 1     2      1345    36    34     7      1422       7.69
#>  2 wt-15-ctl-01 host-t… 2     3      1477    12    40     5      1534      15.4 
#>  3 wt-15-ctl-01 host-t… 3     4      1577    14    12     7      1610      65.7 
#>  4 wt-15-ctl-01 host-t… 4     6      1504     8    15     3      1530      37.6 
#>  5 wt-15-ctl-01 host-t… 6     7      1298    25     0     3      1326     356.  
#>  6 wt-15-ctl-01 host-t… 10    11     1511     7    11     3      1532      58.9 
#>  7 wt-15-ctl-01 host-t… 11    58     1417     2     8     2      1429     177.  
#>  8 wt-15-ctl-01 host-t… 15    16      685    33   290    48      1056       3.44
#>  9 wt-15-ctl-01 host-t… 18    45     1006   728    73   125      1932       2.37
#> 10 wt-15-ctl-01 host-t… 24    25     1852    25    73     8      1958       8.12
#> # ℹ 1,575 more rows
#> # ℹ 7 more variables: log_odds_ratio <dbl>, se_log_or <dbl>, ci_lower <dbl>,
#> #   ci_upper <dbl>, fisher_or <dbl>, p_value <dbl>, p_adjusted <dbl>
```
