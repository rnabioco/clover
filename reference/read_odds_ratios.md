# Read an odds ratios file.

Read a per-position-pair odds ratio file produced by the tRNA sequencing
pipeline. These files contain pairwise modification co-occurrence
statistics.

## Usage

``` r
read_odds_ratios(path)
```

## Arguments

- path:

  Path to a `{sample}.odds_ratios_filtered.tsv.gz` file.

## Value

A tibble with columns including `ref`, `pos1`, `pos2`, `odds_ratio`,
`log_odds_ratio`, `p_value`, and `total_obs`.

## Examples

``` r
path <- clover_example(
  "ecoli/summary/tables/wt-15-ctl-01/wt-15-ctl-01.odds_ratios_filtered.tsv.gz"
)
read_odds_ratios(path)
#> # A tibble: 720 × 16
#>    ref   pos1  pos2    n00   n01   n10   n11 total_obs odds_ratio log_odds_ratio
#>    <chr> <chr> <chr> <dbl> <dbl> <dbl> <dbl>     <dbl>      <dbl>          <dbl>
#>  1 host… 1     2      1345    36    34     7      1422       7.69          2.04 
#>  2 host… 2     3      1477    12    40     5      1534      15.4           2.73 
#>  3 host… 3     4      1577    14    12     7      1610      65.7           4.19 
#>  4 host… 4     6      1504     8    15     3      1530      37.6           3.63 
#>  5 host… 6     7      1298    25     0     3      1326     356.            5.88 
#>  6 host… 10    11     1511     7    11     3      1532      58.9           4.08 
#>  7 host… 11    58     1417     2     8     2      1429     177.            5.18 
#>  8 host… 15    16      685    33   290    48      1056       3.44          1.23 
#>  9 host… 18    45     1006   728    73   125      1932       2.37          0.861
#> 10 host… 24    25     1852    25    73     8      1958       8.12          2.09 
#> # ℹ 710 more rows
#> # ℹ 6 more variables: se_log_or <dbl>, ci_lower <dbl>, ci_upper <dbl>,
#> #   fisher_or <dbl>, p_value <dbl>, p_adjusted <dbl>
```
