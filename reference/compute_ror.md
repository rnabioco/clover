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
#> # A tibble: 8,602 × 6
#>     pos1  pos2 or_numerator or_denominator log_ror     ror
#>    <dbl> <dbl>        <dbl>          <dbl>   <dbl>   <dbl>
#>  1    18    24        -23.0          -23.0    0    1      
#>  2    18    40        -23.0          -23.0    0    1      
#>  3    18    46        -23.0          -23.0    0    1      
#>  4    18    66        -23.0          -23.0    0    1      
#>  5    18    76        -23.0          -23.0    0    1      
#>  6    19    18        -23.0          -23.0    0    1      
#>  7    19    20        -23.0          -17.7   -5.29 0.00506
#>  8    19    21        -23.0          -23.0    0    1      
#>  9    19    22        -23.0          -23.0    0    1      
#> 10    19    24        -23.0          -23.0    0    1      
#> # ℹ 8,592 more rows
```
