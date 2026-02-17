# Compute relative odds ratio at isodecoder level.

Compare isodecoder-level odds ratios between two conditions by computing
the relative odds ratio (ROR) with z-score significance testing.

## Usage

``` r
compute_ror_isodecoder(
  data_num,
  data_den,
  ror_cap = 10,
  p_method = "BH",
  alpha = 0.05
)
```

## Arguments

- data_num:

  A tibble of isodecoder-aggregated odds ratios for the numerator
  condition (from
  [`aggregate_or_isodecoder()`](https://rnabioco.github.io/clover/reference/aggregate_or_isodecoder.md)).

- data_den:

  A tibble of isodecoder-aggregated odds ratios for the denominator
  condition.

- ror_cap:

  Maximum absolute value for the ROR. Values beyond this are capped.
  Default `10`.

- p_method:

  Method for p-value adjustment, passed to
  [`stats::p.adjust()`](https://rdrr.io/r/stats/p.adjust.html). Default
  `"BH"`.

- alpha:

  Significance level for the `significant` column. Default `0.05`.

## Value

A tibble with columns: `isodecoder`, `pos1`, `pos2`, `mean_log_or_num`,
`mean_log_or_den`, `se_num`, `se_den`, `ror`, `ror_se`, `z_score`,
`p_value`, `p_adj`, `ci_lower`, `ci_upper`, and `significant`.

## Examples

``` r
num <- tibble::tibble(
  isodecoder = rep("tRNA-Ala-AGC", 2),
  pos1 = c(20, 34), pos2 = c(34, 58),
  mean_or = c(3.0, 1.5), mean_log_or = c(1.1, 0.4),
  sd_log_or = c(0.2, 0.3), min_pval = c(0.001, 0.01),
  total_reads = c(500, 300), n_copies = c(2, 2)
)
den <- tibble::tibble(
  isodecoder = rep("tRNA-Ala-AGC", 2),
  pos1 = c(20, 34), pos2 = c(34, 58),
  mean_or = c(1.5, 1.2), mean_log_or = c(0.4, 0.18),
  sd_log_or = c(0.15, 0.25), min_pval = c(0.01, 0.05),
  total_reads = c(400, 250), n_copies = c(2, 2)
)
compute_ror_isodecoder(num, den)
#> # A tibble: 2 × 15
#>   isodecoder    pos1  pos2 mean_log_or_num se_num mean_log_or_den se_den   ror
#>   <chr>        <dbl> <dbl>           <dbl>  <dbl>           <dbl>  <dbl> <dbl>
#> 1 tRNA-Ala-AGC    20    34             1.1    0.2            0.4    0.15  0.7 
#> 2 tRNA-Ala-AGC    34    58             0.4    0.3            0.18   0.25  0.22
#> # ℹ 7 more variables: ror_se <dbl>, z_score <dbl>, p_value <dbl>, p_adj <dbl>,
#> #   ci_lower <dbl>, ci_upper <dbl>, significant <lgl>
```
