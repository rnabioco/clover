# Embed a tRNA structure SVG as centered HTML

Reads an SVG file (typically produced by
[`plot_tRNA_structure()`](https://rnabioco.github.io/clover/reference/plot_tRNA_structure.md))
and wraps it in a centering `<div>`, returning an
[`htmltools::HTML()`](https://rstudio.github.io/htmltools/reference/HTML.html)
object suitable for use in R Markdown or Quarto documents.

## Usage

``` r
structure_html(svg_path)
```

## Arguments

- svg_path:

  Path to an SVG file, typically the return value of
  [`plot_tRNA_structure()`](https://rnabioco.github.io/clover/reference/plot_tRNA_structure.md).

## Value

An
[`htmltools::HTML()`](https://rstudio.github.io/htmltools/reference/HTML.html)
object.

## Examples

``` r
if (FALSE) { # \dontrun{
svg <- plot_tRNA_structure("tRNA-Glu-TTC", "Escherichia coli")
structure_html(svg)
} # }
```
