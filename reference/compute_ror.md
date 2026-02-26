# Compute ratio of odds ratios between conditions.

Compare modification co-occurrence between two conditions by computing
the ratio of odds ratios (ROR). Replicates within each condition are
aggregated using the specified function.

## Usage

``` r
compute_ror(
  odds_data,
  condition_col = "condition",
  numerator,
  denominator,
  min_obs = 100,
  agg_fun = mean
)
```

## Arguments

- odds_data:

  A combined tibble of odds ratio data with a `condition` column (or
  column specified by `condition_col`) and `sample_id`.

- condition_col:

  Column name (string) for condition labels. Default `"condition"`.

- numerator:

  Value of `condition_col` for the numerator condition.

- denominator:

  Value of `condition_col` for the denominator condition.

- min_obs:

  Minimum `total_obs` for a pair to be included. Default `100`.

- agg_fun:

  Function to aggregate replicate log odds ratios. Default `mean`.

## Value

A tibble with columns: `pos1`, `pos2`, `or_numerator`, `or_denominator`,
`ror`, and `log_ror`.

## Examples

``` r
results <- read_pipeline_results(
  clover_example("ecoli/config.yaml"),
  types = "odds_ratios"
)
or_data <- results$odds_ratios
or_data$condition <- ifelse(
  grepl("ctl", or_data$sample_id), "ctl", "inf"
)
compute_ror(or_data, numerator = "inf", denominator = "ctl")
#> # A tibble: 525 × 6
#>    pos1  pos2  or_numerator or_denominator log_ror   ror
#>    <chr> <chr>        <dbl>          <dbl>   <dbl> <dbl>
#>  1 1     16          -0.571         -0.899  0.327  1.39 
#>  2 1     2            2.10           2.01   0.0893 1.09 
#>  3 1     20          -0.786         -1.07   0.288  1.33 
#>  4 1     32          -1.34          -1.63   0.285  1.33 
#>  5 1     39          -0.989         -1.08   0.0939 1.10 
#>  6 1     4            1.69           2.13  -0.440  0.644
#>  7 1     40          -1.01          -1.23   0.221  1.25 
#>  8 1     41          -1.04          -1.30   0.260  1.30 
#>  9 1     53          -0.400         -0.833  0.433  1.54 
#> 10 1     55          -0.887         -1.36   0.470  1.60 
#> # ℹ 515 more rows
```
