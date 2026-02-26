# Plot per-tRNA charging ratio differences.

Create a dot plot with error bars showing the difference in charging
ratio between two conditions for each tRNA. Consumes the tibble returned
by
[`compute_charging_diffs()`](https://rnabioco.github.io/clover/reference/compute_charging_diffs.md).

## Usage

``` r
plot_charging_diffs(
  data,
  point_size = 2.5,
  source_col = NULL,
  label_col = "ref",
  shorten = TRUE
)
```

## Arguments

- data:

  A tibble from
  [`compute_charging_diffs()`](https://rnabioco.github.io/clover/reference/compute_charging_diffs.md)
  with at least `ref` (factor), `diff`, and `se_diff` columns.

- point_size:

  Numeric size for
  [`ggplot2::geom_point()`](https://ggplot2.tidyverse.org/reference/geom_point.html).
  Default `2.5`.

- source_col:

  Optional column name (string) for faceting, e.g., `"source"` to
  separate host and phage tRNAs. Default `NULL`.

- label_col:

  Column name (string) to use for y-axis labels. Default `"ref"`.

- shorten:

  Logical; if `TRUE` (default), shorten tRNA names on the y-axis via
  [`shorten_trna_names()`](https://rnabioco.github.io/clover/reference/shorten_trna_names.md).

## Value

A ggplot object.

## Examples

``` r
df <- tibble::tibble(
  ref = forcats::fct_inorder(paste0("tRNA-", 1:5)),
  diff = c(-0.1, -0.05, 0.02, 0.08, 0.15),
  se_diff = rep(0.03, 5)
)
plot_charging_diffs(df)
```
