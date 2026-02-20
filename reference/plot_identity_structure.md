# Plot tRNA structure with identity element overlays

Generates a tRNA cloverleaf structure SVG with aminoacylation identity
elements highlighted as colored outlines. Strong determinants are shown
in red, weak determinants in blue.

## Usage

``` r
plot_identity_structure(
  trna,
  organism,
  sprinzl_coords,
  trna_id = NULL,
  amino_acid = NULL,
  outline_palette = NULL,
  ...
)
```

## Arguments

- trna:

  Character string identifying the tRNA (e.g., `"tRNA-Ala-AGC"`). Use
  [`structure_trnas()`](https://rnabioco.github.io/clover/reference/structure_trnas.md)
  to list available tRNAs.

- organism:

  Character string specifying the organism name (e.g.,
  `"Saccharomyces cerevisiae"`).

- sprinzl_coords:

  A tibble of Sprinzl coordinates as returned by
  [`read_sprinzl_coords()`](https://rnabioco.github.io/clover/reference/read_sprinzl_coords.md).

- trna_id:

  Character string identifying the tRNA in `sprinzl_coords` (e.g.,
  `"nuc-tRNA-Ala-AGC-1-1"`). If `NULL` (default), auto-detected from
  `trna` by converting the anticodon to RNA (T to U) and matching the
  first Sprinzl entry.

- amino_acid:

  Optional 3-letter amino acid code. If `NULL` (default), extracted from
  the `trna` name.

- outline_palette:

  Named character vector of colors keyed by strength (`"strong"`,
  `"weak"`). Default uses red for strong and blue for weak determinants.

- ...:

  Additional arguments passed to
  [`plot_tRNA_structure()`](https://rnabioco.github.io/clover/reference/plot_tRNA_structure.md).

## Value

The path to the annotated SVG file (invisibly).

## Examples

``` r
if (FALSE) { # \dontrun{
coords <- read_sprinzl_coords(
  clover_example("sprinzl/sacCer_global_coords.tsv.gz")
)
plot_identity_structure(
  "tRNA-Ala-AGC", "Saccharomyces cerevisiae", coords
)
} # }
```
