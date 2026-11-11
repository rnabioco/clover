# Plot a query tRNA on a reference cloverleaf, highlighting differences

Draws the bundled cloverleaf for a reference tRNA and overlays how a
query tRNA differs from it, so that several queries share one layout.
Positions are matched through Sprinzl labels (see
[`diff_sprinzl()`](https://rnabioco.github.io/clover/reference/diff_sprinzl.md)):

## Usage

``` r
plot_tRNA_diff(
  trna,
  organism,
  reference_coords,
  reference_id,
  query_coords,
  query_id,
  output = NULL,
  position_markers = TRUE
)
```

## Arguments

- trna, organism:

  The bundled reference structure, as in
  [`plot_tRNA_structure()`](https://rnabioco.github.io/clover/reference/plot_tRNA_structure.md).

- reference_coords, query_coords:

  Tibbles of Sprinzl coordinates as returned by
  [`read_sprinzl_coords()`](https://rnabioco.github.io/clover/reference/read_sprinzl_coords.md).
  The reference positions must correspond to the nucleotides of the
  bundled structure.

- reference_id, query_id:

  The `trna_id` of the reference and query tRNA within their coordinate
  tables.

- output:

  Path for the output SVG file. If `NULL` (default), writes to a
  temporary file.

- position_markers:

  Logical; if `TRUE` (default), draw small grey numbers at Sprinzl
  positions 10, 20, 30, ... Unlike
  [`plot_tRNA_structure()`](https://rnabioco.github.io/clover/reference/plot_tRNA_structure.md),
  these are Sprinzl labels, not sequence indices.

## Value

The path to the annotated SVG file (invisibly).

## Details

- A reference position missing from the query is drawn as an empty grey
  circle with no letter.

- A position with a different base shows the query base, outlined.

- A query base with no reference position cannot sit in the reference
  layout, so it is marked by a small dot on the backbone and listed in
  the legend. Inserted bases at the 5' end are placed at the first
  nucleotide and those at the 3' end at the last.

Nothing in the reference layout moves; bases that match are drawn as in
[`plot_tRNA_structure()`](https://rnabioco.github.io/clover/reference/plot_tRNA_structure.md).

## Examples

``` r
# \donttest{
ecoli <- read_sprinzl_coords(
  clover_example("sprinzl/ecoliK12_global_coords.tsv.gz")
)
t4 <- read_sprinzl_coords(
  clover_example("sprinzl/phageT4_global_coords.tsv.gz")
)
plot_tRNA_diff(
  "tRNA-Gly-TCC",
  "Escherichia coli",
  reference_coords = ecoli,
  reference_id = "tRNA-Gly-UCC-1-1",
  query_coords = t4,
  query_id = "tRNA-Gly-UCC"
)
# }
```
