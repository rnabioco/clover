# Map identity elements to tRNA sequence positions

Converts Sprinzl positions in identity element data to 1-based sequence
positions for a specific tRNA, enabling overlay on structure plots via
the `outlines` parameter of
[`plot_tRNA_structure()`](https://rnabioco.github.io/clover/reference/plot_tRNA_structure.md).

## Usage

``` r
map_identity_to_trna(elements, sprinzl_coords, trna_id)
```

## Arguments

- elements:

  A tibble of identity elements as returned by
  [`identity_elements()`](https://rnabioco.github.io/clover/reference/identity_elements.md).

- sprinzl_coords:

  A tibble of Sprinzl coordinates as returned by
  [`read_sprinzl_coords()`](https://rnabioco.github.io/clover/reference/read_sprinzl_coords.md).

- trna_id:

  Character string identifying the tRNA to map (must match a `trna_id`
  value in `sprinzl_coords`).

## Value

The input `elements` tibble with an added `pos` column containing the
1-based sequence position. Rows where the Sprinzl position does not
exist in the target tRNA (including His G\\\_{-1}\\ at
`sprinzl_pos = -1` and structural determinants with `sprinzl_pos = NA`)
are dropped.

## Examples

``` r
if (FALSE) { # \dontrun{
coords <- read_sprinzl_coords(
  clover_example("sprinzl/sacCer_global_coords.tsv.gz")
)
elems <- identity_elements("Saccharomyces cerevisiae",
  amino_acid = "Ala")
map_identity_to_trna(elems, coords, "nuc-tRNA-Ala-AGC-1-1")
} # }
```
