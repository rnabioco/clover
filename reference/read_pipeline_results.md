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
#> # A tibble: 63,823 × 8
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
#> # ℹ 63,813 more rows
```
