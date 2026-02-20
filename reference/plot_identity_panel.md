# Plot identity elements for multiple tRNAs side by side

Generates a combined SVG showing tRNA cloverleaf structures arranged in
a horizontal row, each annotated with aminoacylation identity elements.
Inspired by Figure 2 of Giege & Eriani (2023).

## Usage

``` r
plot_identity_panel(
  trnas,
  organism,
  sprinzl_coords,
  output = NULL,
  outline_palette = NULL,
  gap = 20,
  ...
)
```

## Arguments

- trnas:

  Character vector of tRNA identifiers (e.g.,
  `c("tRNA-Ala-AGC", "tRNA-Phe-GAA")`).

- organism:

  Character string specifying the organism name.

- sprinzl_coords:

  A tibble of Sprinzl coordinates as returned by
  [`read_sprinzl_coords()`](https://rnabioco.github.io/clover/reference/read_sprinzl_coords.md).

- output:

  Path for the output SVG file. If `NULL` (default), writes to a
  temporary file.

- outline_palette:

  Named character vector of colors keyed by strength. Default uses red
  for strong and blue for weak.

- gap:

  Horizontal gap in SVG units between panels. Default 20.

- ...:

  Additional arguments passed to
  [`plot_identity_structure()`](https://rnabioco.github.io/clover/reference/plot_identity_structure.md).

## Value

The path to the combined SVG file (invisibly).

## Examples

``` r
if (FALSE) { # \dontrun{
coords <- read_sprinzl_coords(
  clover_example("sprinzl/sacCer_global_coords.tsv.gz")
)
plot_identity_panel(
  c("tRNA-Ala-AGC", "tRNA-Asp-GTC",
    "tRNA-Phe-GAA", "tRNA-His-GTG"),
  "Saccharomyces cerevisiae",
  coords
)
} # }
```
