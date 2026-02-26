# Plot per-tRNA charging ratios.

Create a box-and-jitter plot of raw charging ratios grouped by
condition. Expects the tibble stored in `metadata(se)$charging_ratios`.

## Usage

``` r
plot_charging_ratios(
  data,
  group_col = "condition",
  facet_col = NULL,
  shorten = TRUE,
  point_size = 0.8,
  point_alpha = 0.4
)
```

## Arguments

- data:

  A tibble with at least `ref`, `charging_ratio`, and the column named
  by `group_col`.

- group_col:

  Column name (string) for the x-axis grouping variable (e.g.,
  `"condition"`). Default `"condition"`.

- facet_col:

  Optional column name (string) for faceting (e.g., `"strain"`). Default
  `NULL`.

- shorten:

  Logical; if `TRUE` (default), shorten tRNA names on the y-axis via
  [`shorten_trna_names()`](https://rnabioco.github.io/clover/reference/shorten_trna_names.md).

- point_size:

  Numeric size for
  [`ggplot2::geom_jitter()`](https://ggplot2.tidyverse.org/reference/geom_jitter.html).
  Default `0.8`.

- point_alpha:

  Numeric alpha for jittered points. Default `0.4`.

## Value

A ggplot object.

## Examples

``` r
df <- tibble::tibble(
  ref = rep(paste0("tRNA-Ala-AGC-", 1:3, "-1"), each = 6),
  condition = rep(c("ctl", "inf"), each = 3, times = 3),
  charging_ratio = runif(18, 0.3, 0.9)
)
plot_charging_ratios(df)
```
