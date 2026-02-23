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
coords <- read_sprinzl_coords(
  clover_example("sprinzl/sacCer_global_coords.tsv.gz")
)
elems <- identity_elements(
  "Saccharomyces cerevisiae",
  amino_acid = "Ala"
)
map_identity_to_trna(elems, coords, "nuc-tRNA-Ala-AGC-1-1")
#> # A tibble: 5 × 13
#>   amino_acid domain  sprinzl_pos nucleotide region       type  strength pair_pos
#>   <chr>      <chr>         <int> <chr>      <chr>        <chr> <chr>       <int>
#> 1 Ala        Eukarya          73 A          discriminat… dete… strong         NA
#> 2 Ala        Eukarya           3 G          acceptor_st… dete… strong         70
#> 3 Ala        Eukarya          70 U          acceptor_st… dete… strong          3
#> 4 Ala        Eukarya           3 G          acceptor_st… anti… NA             70
#> 5 Ala        Eukarya          70 U          acceptor_st… anti… NA              3
#> # ℹ 5 more variables: pair_type <chr>, against_aars <chr>, universal <lgl>,
#> #   description <chr>, pos <dbl>
```
