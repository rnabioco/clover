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
if (FALSE) { # \dontrun{
paths <- c(
  ctl1 = "ctl1.charging.cpm.tsv.gz",
  ctl2 = "ctl2.charging.cpm.tsv.gz",
  inf1 = "inf1.charging.cpm.tsv.gz",
  inf2 = "inf2.charging.cpm.tsv.gz"
)
charging <- read_charging_multi(paths)
charging$condition <- ifelse(
  grepl("ctl", charging$sample_id), "ctl", "inf"
)
diffs <- compute_charging_diffs(
  charging,
  numerator = "inf",
  denominator = "ctl"
)
} # }
```
