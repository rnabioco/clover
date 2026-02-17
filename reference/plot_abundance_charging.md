# Plot abundance changes versus charging ratio changes.

Creates a scatter plot comparing tRNA abundance changes (from DESeq2)
with charging ratio changes on a single plot. Significant points are
colored by quadrant and labeled with
[`ggrepel::geom_text_repel()`](https://ggrepel.slowkow.com/reference/geom_text_repel.html).

## Usage

``` r
plot_abundance_charging(
  deseq_res,
  charging_diffs,
  lab_col = "tRNA",
  padj_cutoff = 0.05,
  max_overlaps = 20,
  point_size = 2,
  label_size = 3
)
```

## Arguments

- deseq_res:

  A tibble from
  [`tidy_deseq_results()`](https://rnabioco.github.io/clover/reference/tidy_deseq_results.md)
  with at least `tRNA`, `log2FoldChange`, and `padj` columns.

- charging_diffs:

  A tibble from
  [`compute_charging_diffs()`](https://rnabioco.github.io/clover/reference/compute_charging_diffs.md)
  with at least `tRNA` and `diff` columns.

- lab_col:

  Column name (string) used for point labels. Default `"tRNA"`.

- padj_cutoff:

  Numeric; significance threshold for `padj`. Default `0.05`.

- max_overlaps:

  Maximum number of overlapping labels passed to
  [`ggrepel::geom_text_repel()`](https://ggrepel.slowkow.com/reference/geom_text_repel.html).
  Default `20`.

- point_size:

  Numeric size for
  [`ggplot2::geom_point()`](https://ggplot2.tidyverse.org/reference/geom_point.html).
  Default `2`.

- label_size:

  Numeric size for
  [`ggrepel::geom_text_repel()`](https://ggrepel.slowkow.com/reference/geom_text_repel.html).
  Default `3`.

## Value

A ggplot object.

## Examples

``` r
deseq_res <- tibble::tibble(
  tRNA = paste0("tRNA-", 1:6),
  log2FoldChange = c(1, -1, 0.5, -0.5, 2, -2),
  padj = c(0.01, 0.02, 0.5, 0.6, 0.001, 0.003)
)
charging_diffs <- tibble::tibble(
  tRNA = paste0("tRNA-", 1:6),
  diff = c(0.1, -0.1, 0.05, -0.05, -0.2, 0.15),
  se_diff = rep(0.03, 6)
)
plot_abundance_charging(deseq_res, charging_diffs)
```
