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
bcerror_rds <- readRDS(clover_example("ecoli/bcerror_summary.rds"))
bcerror_delta <- compute_bcerror_delta(
  bcerror_rds$bcerror_summary, delta = wt - tb
)
sprinzl <- read_sprinzl_coords(
  clover_example("sprinzl/ecoliK12_global_coords.tsv.gz")
)
prep_mod_heatmap(bcerror_delta, sprinzl_coords = sprinzl)
#> # A tibble: 776 × 9
#>    ref             pos     tb     wt    delta trna_id sprinzl_label global_index
#>    <chr>         <dbl>  <dbl>  <dbl>    <dbl> <chr>   <fct>                <dbl>
#>  1 tRNA-Asn-GTT…     1 0.0334 0.0405  7.11e-3 tRNA-A… 1                        1
#>  2 tRNA-Asn-GTT…     2 0.0660 0.0896  2.36e-2 tRNA-A… 2                        3
#>  3 tRNA-Asn-GTT…     3 0.0587 0.0809  2.21e-2 tRNA-A… 3                        4
#>  4 tRNA-Asn-GTT…     4 0.195  0.201   5.84e-3 tRNA-A… 4                        8
#>  5 tRNA-Asn-GTT…     5 0.258  0.288   2.92e-2 tRNA-A… 5                        9
#>  6 tRNA-Asn-GTT…     6 0.117  0.134   1.70e-2 tRNA-A… 6                       11
#>  7 tRNA-Asn-GTT…     7 0.125  0.133   7.20e-3 tRNA-A… 7                       12
#>  8 tRNA-Asn-GTT…     8 0.0644 0.0706  6.17e-3 tRNA-A… 8                       13
#>  9 tRNA-Asn-GTT…     9 0.0623 0.0619 -3.87e-4 tRNA-A… 9                       14
#> 10 tRNA-Asn-GTT…    10 0.0295 0.0266 -2.95e-3 tRNA-A… 10                      15
#> # ℹ 766 more rows
#> # ℹ 1 more variable: trna_label <chr>
```
