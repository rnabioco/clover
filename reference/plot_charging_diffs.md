# Plot per-tRNA charging ratio differences.

Create a dot plot with error bars showing the difference in charging
ratio between two conditions for each tRNA. Consumes the tibble returned
by
[`compute_charging_diffs()`](https://rnabioco.github.io/clover/reference/compute_charging_diffs.md).

## Usage

``` r
plot_charging_diffs(data, point_size = 2.5)
```

## Arguments

- data:

  A tibble from
  [`compute_charging_diffs()`](https://rnabioco.github.io/clover/reference/compute_charging_diffs.md)
  with at least `tRNA` (factor), `diff`, and `se_diff` columns.

- point_size:

  Numeric size for
  [`ggplot2::geom_point()`](https://ggplot2.tidyverse.org/reference/geom_point.html).
  Default `2.5`.

## Value

A ggplot object.

## Examples

``` r
df <- tibble::tibble(
  tRNA = forcats::fct_inorder(paste0("tRNA-", 1:5)),
  diff = c(-0.1, -0.05, 0.02, 0.08, 0.15),
  se_diff = rep(0.03, 5)
)
plot_charging_diffs(df)
```
