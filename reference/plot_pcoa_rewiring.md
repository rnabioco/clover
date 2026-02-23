# Plot PCoA of tRNA rewiring scores.

Create a scatter plot of PCoA coordinates from
[`perform_pcoa()`](https://rnabioco.github.io/clover/reference/perform_pcoa.md),
with points sized by the number of non-zero comparisons and colored by
Euclidean rewiring magnitude. The top `n_label` isodecoders are labeled
with
[`ggrepel::geom_text_repel()`](https://ggrepel.slowkow.com/reference/geom_text_repel.html).

## Usage

``` r
plot_pcoa_rewiring(
  pcoa_result,
  rewiring_scores,
  title = "PCoA of tRNA rewiring",
  n_label = 10
)
```

## Arguments

- pcoa_result:

  A list from
  [`perform_pcoa()`](https://rnabioco.github.io/clover/reference/perform_pcoa.md)
  with elements `coordinates` and `variance_explained`.

- rewiring_scores:

  A tibble from
  [`calculate_rewiring_scores()`](https://rnabioco.github.io/clover/reference/calculate_rewiring_scores.md)
  with at least `isodecoder`, `euclidean_magnitude`, and `n_nonzero`
  columns.

- title:

  Plot title. Default `"PCoA of tRNA rewiring"`.

- n_label:

  Number of top isodecoders to label. Default `10`.

## Value

A ggplot object.

## Examples

``` r
# \donttest{
mat <- matrix(
  c(1.5, -0.8, 0.3, 2.1, 0.5, -1.2),
  nrow = 3,
  dimnames = list(
    c("tRNA-Ala", "tRNA-Gly", "tRNA-Ser"),
    c("20_vs_34", "34_vs_58")
  )
)
scores <- calculate_rewiring_scores(mat)
pcoa <- perform_pcoa(mat)
plot_pcoa_rewiring(pcoa, scores)

# }
```
