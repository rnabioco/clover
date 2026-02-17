# Plot a delta-signal modification heatmap.

Create a diverging heatmap of modification signal changes (e.g., mutant
minus wild-type) across tRNA families and Sprinzl positions. Rows can be
optionally clustered using Ward's D2 hierarchical clustering.

## Usage

``` r
plot_mod_heatmap(
  data,
  value_col = "value",
  ref_col = "ref",
  cluster = TRUE,
  color_limits = c(-0.25, 0.25),
  color_low = "#0072B2",
  color_high = "#D55E00",
  na_value = "gray80",
  square = TRUE
)
```

## Arguments

- data:

  A data frame with at least three columns: one for tRNA
  family/reference (y-axis), one for Sprinzl position labels (x-axis),
  and one for the fill value.

- value_col:

  Column name (string) for fill values. Default `"value"`.

- ref_col:

  Column name (string) for tRNA families (y-axis). Default `"ref"`.

- cluster:

  Logical; cluster rows with Ward's D2? Default `TRUE`.

- color_limits:

  Numeric vector of length 2 giving symmetric limits for the color
  scale. Default `c(-0.25, 0.25)`.

- color_low:

  Color for negative values. Default `"#0072B2"` (blue).

- color_high:

  Color for positive values. Default `"#D55E00"` (red).

- na_value:

  Color for missing positions. Default `"gray80"`.

- square:

  Logical; use `coord_fixed(ratio = 1)`? Default `TRUE`.

## Value

A ggplot object.

## Examples

``` r
df <- tidyr::expand_grid(
  ref = paste0("tRNA-", c("Ala", "Gly", "Ser")),
  sprinzl_label = as.character(1:10)
)
df$value <- rnorm(nrow(df), sd = 0.1)
plot_mod_heatmap(df)
```
