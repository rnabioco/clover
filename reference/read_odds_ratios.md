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

  Path to a `{sample}.odds_ratios.tsv.gz` file.

## Value

A tibble with columns including `ref`, `pos1`, `pos2`, `odds_ratio`,
`log_odds_ratio`, `p_value`, and `total_obs`.

## Examples

``` r
path <- clover_example(
  "ecoli/summary/tables/wt-15-ctl-01/wt-15-ctl-01.odds_ratios.tsv.gz"
)
read_odds_ratios(path)
#> # A tibble: 11,310 × 7
#>    ref                    pos1  pos2 odds_ratio log_odds_ratio p_value total_obs
#>    <chr>                 <dbl> <dbl>      <dbl>          <dbl>   <dbl>     <dbl>
#>  1 host-tRNA-Asp-GTC-1-1    20    23          0          -23.0       1      1488
#>  2 host-tRNA-Asp-GTC-1-1    20    26          0          -23.0       1      1488
#>  3 host-tRNA-Asp-GTC-1-1    20    28          0          -23.0       1      1488
#>  4 host-tRNA-Asp-GTC-1-1    20    31          0          -23.0       1      1488
#>  5 host-tRNA-Asp-GTC-1-1    20    32          0          -23.0       1      1488
#>  6 host-tRNA-Asp-GTC-1-1    20    35          0          -23.0       1      1488
#>  7 host-tRNA-Asp-GTC-1-1    20    36          0          -23.0       1      1488
#>  8 host-tRNA-Asp-GTC-1-1    20    37          0          -23.0       1      1488
#>  9 host-tRNA-Asp-GTC-1-1    20    39          0          -23.0       1      1488
#> 10 host-tRNA-Asp-GTC-1-1    20    40          0          -23.0       1      1488
#> # ℹ 11,300 more rows
```
