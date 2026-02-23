# Plot tRNA secondary structure with modifications and linkages

Reads a bundled tRNA cloverleaf SVG and overlays modification
highlights, outline circles, and circuit linkage arcs. Modifications are
shown as colored filled circles behind nucleotide letters; outlines are
shown as colored circle borders; linkages are drawn as Bezier curve arcs
between position pairs.

## Usage

``` r
plot_tRNA_structure(
  trna,
  organism,
  modifications = NULL,
  outlines = NULL,
  linkages = NULL,
  output = NULL,
  mod_palette = NULL,
  outline_palette = NULL,
  text_colors = NULL,
  position_markers = TRUE,
  linkage_palette = c("#0072B2", "#D55E00")
)
```

## Arguments

- trna:

  Character string identifying the tRNA (e.g., `"tRNA-Ala-GGC"`). Use
  [`structure_trnas()`](https://rnabioco.github.io/clover/reference/structure_trnas.md)
  to list available tRNAs.

- organism:

  Character string specifying the organism name (e.g.,
  `"Escherichia coli"`).

- modifications:

  A tibble with columns `pos` (1-based position in the tRNA sequence)
  and `mod1` (short modification name, e.g., `"m1A"`). Output of
  [`modomics_mods()`](https://rnabioco.github.io/clover/reference/modomics_mods.md)
  works directly after filtering to the tRNA of interest.

- outlines:

  A tibble with columns `pos` (1-based position) and `group` (category
  name for palette lookup). Draws circle outlines (stroke only, no fill)
  around each nucleotide.

- linkages:

  A tibble with columns `pos1`, `pos2`, and optionally `value` (e.g.,
  log odds ratio) for coloring arcs. If a `log_odds_ratio` column is
  present and `value` is not, it is automatically used as `value`, so
  output of
  [`clean_odds_ratios()`](https://rnabioco.github.io/clover/reference/clean_odds_ratios.md)
  or
  [`filter_linkages()`](https://rnabioco.github.io/clover/reference/filter_linkages.md)
  works directly.

- output:

  Path for the output SVG file. If `NULL` (default), writes to a
  temporary file.

- mod_palette:

  Named character vector of colors keyed by modification short name. If
  `NULL`, uses a default palette.

- outline_palette:

  Named character vector of colors keyed by outline group name. If
  `NULL`, uses `"#333333"` for all.

- text_colors:

  A tibble with columns `pos` (1-based position) and `color` (hex color
  string). Changes the nucleotide letter color at specified positions.
  Unspecified positions keep the default color.

- position_markers:

  Logical; if `TRUE` (default), draw small grey position numbers every
  10 nucleotides around the cloverleaf to help orient readers.

- linkage_palette:

  Character vector of length 2 giving the colors for negative
  (exclusive) and positive (co-occurring) linkage values. Default
  `c("#0072B2", "#D55E00")` (blue for exclusive, vermillion for
  co-occurring). Stroke width encodes the magnitude of the value.

## Value

The path to the annotated SVG file (invisibly).

## Examples

``` r
# \donttest{
plot_tRNA_structure("tRNA-Glu-TTC", "Escherichia coli")
# }
```
