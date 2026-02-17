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
  min_obs = 50,
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

  Minimum `total_obs` for a pair to be included. Default `50`.

- agg_fun:

  Function to aggregate replicate log odds ratios. Default `mean`.

## Value

A tibble with columns: `pos1`, `pos2`, `or_numerator`, `or_denominator`,
`ror`, and `log_ror`.

## Examples

``` r
if (FALSE) { # \dontrun{
combined_or <- read_odds_ratios_multi(paths)
combined_or$condition <- ifelse(
  grepl("wt", combined_or$sample_id), "wt", "mut"
)
ror <- compute_ror(combined_or, numerator = "mut", denominator = "wt")
} # }
```
