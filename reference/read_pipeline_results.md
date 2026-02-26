# Read pipeline results for all samples.

Convenience function that reads a pipeline config, discovers output
files, and loads data for the requested result types. Each result type
is returned as a combined tibble with a `sample_id` column.

## Usage

``` r
read_pipeline_results(
  config_path,
  types = c("charging", "bcerror", "odds_ratios")
)
```

## Arguments

- config_path:

  Path to a Snakemake `config.yaml` file.

- types:

  Character vector of result types to load. Valid values: `"charging"`,
  `"bcerror"`, `"odds_ratios"`.

## Value

A named list of tibbles, one per requested type. Each tibble includes a
`sample_id` column identifying the source sample.

## Examples

``` r
results <- read_pipeline_results(clover_example("ecoli/config.yaml"))
results$charging
#> # A tibble: 1,140 × 6
#>    sample_id    ref    counts_charged counts_uncharged cpm_charged cpm_uncharged
#>    <chr>        <chr>           <dbl>            <dbl>       <dbl>         <dbl>
#>  1 wt-15-ctl-01 host-…            174              124        270.          193.
#>  2 wt-15-ctl-01 host-…           1173             3220       1824.         5006.
#>  3 wt-15-ctl-01 host-…            818              807       1272.         1255.
#>  4 wt-15-ctl-01 host-…           1205             3186       1873.         4953.
#>  5 wt-15-ctl-01 host-…            351              167        546.          260.
#>  6 wt-15-ctl-01 host-…           1390             2395       2161.         3723.
#>  7 wt-15-ctl-01 host-…            615              214        956.          333.
#>  8 wt-15-ctl-01 host-…           4739            13203       7367.        20525.
#>  9 wt-15-ctl-01 host-…           1167             1162       1814.         1806.
#> 10 wt-15-ctl-01 host-…           1321             2330       2054.         3622.
#> # ℹ 1,130 more rows
results$odds_ratios
#> # A tibble: 5,432 × 17
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
#> # ℹ 5,422 more rows
#> # ℹ 7 more variables: log_odds_ratio <dbl>, se_log_or <dbl>, ci_lower <dbl>,
#> #   ci_upper <dbl>, fisher_or <dbl>, p_value <dbl>, p_adjusted <dbl>
```
