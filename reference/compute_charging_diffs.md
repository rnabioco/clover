# Compute charging ratio differences between conditions.

Compare per-tRNA charging ratios (charged / total) between two
conditions, summarizing replicates with mean and standard error and
computing the between-condition difference with propagated SE.

## Usage

``` r
compute_charging_diffs(
  charging_data,
  condition_col = "condition",
  numerator,
  denominator,
  min_count = 50,
  n_top = NULL
)
```

## Arguments

- charging_data:

  A tibble from
  [`read_charging_multi()`](https://rnabioco.github.io/clover/reference/read_charging_multi.md)
  with an added condition column. Must contain `ref`, `counts_charged`,
  `counts_uncharged`, `sample_id`, and the column named by
  `condition_col`.

- condition_col:

  Column name (string) for condition labels. Default `"condition"`.

- numerator:

  Value of `condition_col` for the numerator condition.

- denominator:

  Value of `condition_col` for the denominator condition.

- min_count:

  Minimum total reads (charged + uncharged) per tRNA per sample to
  include. Default `50`.

- n_top:

  If non-`NULL`, keep only the top `n_top` tRNAs by total abundance
  across all samples. Default `NULL` (keep all).

## Value

A tibble with columns:

- ref:

  tRNA identifier (factor ordered by `diff`).

- ratio_numerator:

  Mean charging ratio for the numerator condition.

- ratio_denominator:

  Mean charging ratio for the denominator condition.

- se_numerator:

  Standard error of the numerator ratio.

- se_denominator:

  Standard error of the denominator ratio.

- diff:

  Difference: `ratio_numerator - ratio_denominator`.

- se_diff:

  Propagated SE: `sqrt(se_numerator^2 + se_denominator^2)`.

## Examples

``` r
config_path <- clover_example("ecoli/config.yaml")
results <- read_pipeline_results(config_path, types = "charging")
charging <- results$charging
charging$condition <- ifelse(
  grepl("ctl", charging$sample_id), "ctl", "inf"
)
compute_charging_diffs(
  charging,
  numerator = "inf",
  denominator = "ctl"
)
#> # A tibble: 88 × 7
#>    ref     ratio_numerator ratio_denominator se_numerator se_denominator    diff
#>    <fct>             <dbl>             <dbl>        <dbl>          <dbl>   <dbl>
#>  1 host-t…           0.552             0.581       0.0503         0.0177 -0.0282
#>  2 host-t…           0.462             0.509       0.0550         0.0109 -0.0477
#>  3 host-t…           0.533             0.601       0.0572         0.0383 -0.0679
#>  4 host-t…           0.622             0.659       0.0359         0.0547 -0.0362
#>  5 host-t…           0.447             0.521       0.0546         0.0101 -0.0748
#>  6 host-t…           0.585             0.679       0.0641         0.0292 -0.0940
#>  7 host-t…           0.624             0.706       0.0530         0.0465 -0.0826
#>  8 host-t…           0.601             0.686       0.0541         0.0313 -0.0849
#>  9 host-t…           0.594             0.705       0.0510         0.0173 -0.111 
#> 10 host-t…           0.330             0.436       0.0491         0.0419 -0.106 
#> # ℹ 78 more rows
#> # ℹ 1 more variable: se_diff <dbl>
```
