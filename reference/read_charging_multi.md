# Read charging CPM files for multiple samples.

Read and combine charging CPM files from multiple samples into a single
tibble with a `sample_id` column.

## Usage

``` r
read_charging_multi(paths)
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
    "ecoli/summary/tables/wt-15-ctl-01/wt-15-ctl-01.charging.cpm.tsv.gz"
  ),
  "wt-15-ctl-02" = clover_example(
    "ecoli/summary/tables/wt-15-ctl-02/wt-15-ctl-02.charging.cpm.tsv.gz"
  )
)
read_charging_multi(paths)
#> # A tibble: 380 × 6
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
#> # ℹ 370 more rows
```
