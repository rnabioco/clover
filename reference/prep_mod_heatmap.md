# Prepare data for [`plot_mod_heatmap()`](https://rnabioco.github.io/clover/reference/plot_mod_heatmap.md).

Join Sprinzl coordinates to a bcerror delta tibble, optionally annotate
known modifications, order positions by canonical Sprinzl index, and
shorten tRNA names for display.

## Usage

``` r
prep_mod_heatmap(
  data,
  value_col = "delta",
  ref_col = "ref",
  sprinzl_coords,
  mods = NULL,
  strip_prefix = "^host-",
  shorten_labels = TRUE
)
```

## Arguments

- data:

  A tibble with at least `ref`, `pos`, and a value column (e.g., `delta`
  from
  [`compute_bcerror_delta()`](https://rnabioco.github.io/clover/reference/compute_bcerror_delta.md)).

- value_col:

  Column name (string) for the fill value. Default `"delta"`.

- ref_col:

  Column name (string) for the tRNA reference. Default `"ref"`.

- sprinzl_coords:

  A tibble of Sprinzl coordinates from
  [`read_sprinzl_coords()`](https://rnabioco.github.io/clover/reference/read_sprinzl_coords.md),
  with at least `trna_id`, `pos`, `sprinzl_label`, and `global_index`
  columns.

- mods:

  Optional tibble of modification annotations (e.g., from
  [`modomics_mods()`](https://rnabioco.github.io/clover/reference/modomics_mods.md))
  with at least `ref` and `pos` columns. When provided, a logical
  `has_mod` column is added.

- strip_prefix:

  Regex pattern to strip from `ref` before matching Sprinzl coordinates.
  Default `"^host-"`.

- shorten_labels:

  Logical; if `TRUE` (default), create a `trna_label` column with
  shortened names (e.g., `"tRNA-Glu-TTC-1-1"` becomes `"Glu-TTC"`).

## Value

A tibble ready for
[`plot_mod_heatmap()`](https://rnabioco.github.io/clover/reference/plot_mod_heatmap.md),
with `sprinzl_label` as an ordered factor and optionally `has_mod` and
`trna_label` columns.

## Examples

``` r
if (FALSE) { # \dontrun{
bcerror_delta <- compute_bcerror_delta(bcerror_summary, delta = wt - tb)
sprinzl <- read_sprinzl_coords(
  clover_example("sprinzl/ecoliK12_global_coords.tsv.gz")
)
mods <- modomics_mods(trna_fasta, organism = "Escherichia coli")
heatmap_data <- prep_mod_heatmap(
  bcerror_delta,
  sprinzl_coords = sprinzl,
  mods = mods
)
} # }
```
