# Convert a tRNA structure SVG to PNG

Renders an SVG file (typically produced by
[`plot_tRNA_structure()`](https://rnabioco.github.io/clover/reference/plot_tRNA_structure.md))
to a PNG bitmap. Requires the
[rsvg](https://cran.r-project.org/package=rsvg) package.

## Usage

``` r
structure_to_png(svg_path, output = NULL, width = NULL, height = NULL)
```

## Arguments

- svg_path:

  Path to an SVG file, typically the return value of
  [`plot_tRNA_structure()`](https://rnabioco.github.io/clover/reference/plot_tRNA_structure.md).

- output:

  Path for the output PNG file. If `NULL` (default), replaces the `.svg`
  extension with `.png`.

- width:

  Width of the output PNG in pixels. If `NULL` (default), uses the
  intrinsic SVG width.

- height:

  Height of the output PNG in pixels. If `NULL` (default), uses the
  intrinsic SVG height.

## Value

The path to the PNG file (invisibly).

## Examples

``` r
# \donttest{
svg <- plot_tRNA_structure("tRNA-Glu-TTC", "Escherichia coli")
structure_to_png(svg)
# }
```
