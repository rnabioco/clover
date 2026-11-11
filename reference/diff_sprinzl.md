# Compare a query tRNA to a reference by Sprinzl label.

Aligns two tRNAs through their Sprinzl position labels and classifies
each reference position as `"match"`, `"substitution"` (same label,
different base) or `"deletion"` (label absent from the query). Query
positions whose label is absent from the reference are returned
separately as insertions.

## Usage

``` r
diff_sprinzl(reference_coords, reference_id, query_coords, query_id)
```

## Arguments

- reference_coords, query_coords:

  Tibbles of Sprinzl coordinates as returned by
  [`read_sprinzl_coords()`](https://rnabioco.github.io/clover/reference/read_sprinzl_coords.md).

- reference_id, query_id:

  The `trna_id` of the reference and query tRNA within their coordinate
  tables.

## Value

A list with two tibbles:

- `reference`: one row per reference position with columns `pos`,
  `sprinzl_label`, `ref_base`, `query_pos`, `query_base` and `status`.

- `insertions`: one row per query-only position with columns
  `query_pos`, `sprinzl_label`, `query_base`, `after_pos` and
  `after_label`, the reference position (sequence index and Sprinzl
  label) the inserted base follows. Both are `0`/`NA` at the 5' end.

## Examples

``` r
ecoli <- read_sprinzl_coords(
  clover_example("sprinzl/ecoliK12_global_coords.tsv.gz")
)
t4 <- read_sprinzl_coords(
  clover_example("sprinzl/phageT4_global_coords.tsv.gz")
)
d <- diff_sprinzl(ecoli, "tRNA-Gly-UCC-1-1", t4, "tRNA-Gly-UCC")
d$reference[d$reference$status != "match", ]
#> # A tibble: 13 × 6
#>      pos sprinzl_label ref_base query_pos query_base status      
#>    <int> <chr>         <chr>        <int> <chr>      <chr>       
#>  1     5 5             G                5 A          substitution
#>  2     6 6             C                6 U          substitution
#>  3    19 20            C               19 U          substitution
#>  4    20 20a           U               NA NA         deletion    
#>  5    31 31            C               30 A          substitution
#>  6    39 39            G               38 U          substitution
#>  7    49 50            C               48 U          substitution
#>  8    51 52            G               50 A          substitution
#>  9    61 62            C               60 U          substitution
#> 10    63 64            G               62 A          substitution
#> 11    64 65            C               63 U          substitution
#> 12    66 67            G               65 A          substitution
#> 13    67 68            C               66 U          substitution
```
