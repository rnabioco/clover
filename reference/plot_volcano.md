# Plot a volcano plot of differential expression results.

Creates a volcano plot from the tibble returned by
[`tidy_deseq_results()`](https://rnabioco.github.io/clover/reference/tidy_deseq_results.md).
Significant points are labeled with
[`ggrepel::geom_text_repel()`](https://ggrepel.slowkow.com/reference/geom_text_repel.html).

## Usage

``` r
plot_volcano(
  data,
  lab_col = "tRNA",
  padj_cutoff = 0.05,
  max_overlaps = 20,
  point_size = 1.5,
  label_size = 3,
  sig_color = "#D55E00",
  nonsig_color = "grey60"
)
```

## Arguments

- data:

  A tibble from
  [`tidy_deseq_results()`](https://rnabioco.github.io/clover/reference/tidy_deseq_results.md)
  with at least `log2FoldChange`, `pvalue`, and `significant` columns.

- lab_col:

  Column name (string) used for point labels. Default `"tRNA"`.

- padj_cutoff:

  Numeric; draws a dashed horizontal line at `-log10(padj_cutoff)`.
  Default `0.05`.

- max_overlaps:

  Maximum number of overlapping labels passed to
  [`ggrepel::geom_text_repel()`](https://ggrepel.slowkow.com/reference/geom_text_repel.html).
  Default `20`.

- point_size:

  Numeric size for
  [`ggplot2::geom_point()`](https://ggplot2.tidyverse.org/reference/geom_point.html).
  Default `1.5`.

- label_size:

  Numeric size for
  [`ggrepel::geom_text_repel()`](https://ggrepel.slowkow.com/reference/geom_text_repel.html).
  Default `3`.

- sig_color:

  Color for significant points. Default `"#D55E00"`.

- nonsig_color:

  Color for non-significant points. Default `"grey60"`.

## Value

A ggplot object.

## Examples

``` r
res <- tibble::tibble(
  tRNA = paste0("tRNA-", 1:10),
  log2FoldChange = rnorm(10),
  pvalue = c(rep(0.001, 3), rep(0.5, 7)),
  padj = c(rep(0.01, 3), rep(0.8, 7)),
  significant = c(rep(TRUE, 3), rep(FALSE, 7))
)
plot_volcano(res)
```
