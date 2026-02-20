# Retrieve tRNA aminoacylation identity elements

Returns experimentally validated identity elements (determinants and
antideterminants) for tRNA aminoacylation, based on data from Giege &
Eriani (2023). Elements are defined per amino acid family using Sprinzl
numbering.

## Usage

``` r
identity_elements(
  organism,
  amino_acid = NULL,
  type = c("both", "determinant", "antideterminant")
)
```

## Arguments

- organism:

  Character string specifying the organism name (e.g.,
  `"Escherichia coli"`, `"Saccharomyces cerevisiae"`, `"Homo sapiens"`).
  Use
  [`identity_organisms()`](https://rnabioco.github.io/clover/reference/identity_organisms.md)
  to list supported organisms.

- amino_acid:

  Optional character vector of 3-letter amino acid codes to filter by
  (e.g., `"Ala"`, `"Arg"`). If `NULL` (default), returns elements for
  all amino acids.

- type:

  Which element types to return: `"both"` (default), `"determinant"`, or
  `"antideterminant"`.

## Value

A tibble with columns:

- `amino_acid`: 3-letter amino acid code

- `domain`: evolutionary domain (`"Bacteria"` or `"Eukarya"`)

- `sprinzl_pos`: Sprinzl position number (integer; `NA` for structural
  determinants like the long variable arm; `-1` for the His G\\\_{-1}\\
  position)

- `nucleotide`: expected base (A/C/G/U, or `NA` if variable)

- `region`: tRNA structural region name

- `type`: `"determinant"` or `"antideterminant"`

- `strength`: `"strong"` or `"weak"` (determinants only)

- `pair_pos`: partner Sprinzl position for base pairs (`NA` for
  single-nucleotide elements)

- `pair_type`: `"wc"` or `"wobble"` (`NA` for unpaired)

- `against_aars`: aaRS blocked (antideterminants only)

- `universal`: `TRUE` if conserved across all domains

- `description`: human-readable description

## References

Giege R, Eriani G (2023). "The tRNA identity landscape for
aminoacylation and beyond." *Nucleic Acids Research*, 51(4), 1528–1570.
[doi:10.1093/nar/gkad007](https://doi.org/10.1093/nar/gkad007)

## Examples

``` r
# All identity elements for E. coli
identity_elements("Escherichia coli")
#> # A tibble: 135 × 12
#>    amino_acid domain   sprinzl_pos nucleotide region     type  strength pair_pos
#>    <chr>      <chr>          <int> <chr>      <chr>      <chr> <chr>       <int>
#>  1 Ala        Bacteria          73 A          discrimin… dete… strong         NA
#>  2 Ala        Bacteria           3 G          acceptor_… dete… strong         70
#>  3 Ala        Bacteria          70 U          acceptor_… dete… strong          3
#>  4 Ala        Bacteria           4 G          acceptor_… dete… weak           69
#>  5 Ala        Bacteria          69 C          acceptor_… dete… weak            4
#>  6 Ala        Bacteria          20 G          d_arm      dete… weak           NA
#>  7 Arg        Bacteria          73 NA         discrimin… dete… strong         NA
#>  8 Arg        Bacteria           2 G          acceptor_… dete… weak           71
#>  9 Arg        Bacteria          71 C          acceptor_… dete… weak            2
#> 10 Arg        Bacteria           3 C          acceptor_… dete… weak           70
#> # ℹ 125 more rows
#> # ℹ 4 more variables: pair_type <chr>, against_aars <chr>, universal <lgl>,
#> #   description <chr>

# Ala determinants only
identity_elements("Escherichia coli", amino_acid = "Ala",
  type = "determinant")
#> # A tibble: 6 × 12
#>   amino_acid domain   sprinzl_pos nucleotide region      type  strength pair_pos
#>   <chr>      <chr>          <int> <chr>      <chr>       <chr> <chr>       <int>
#> 1 Ala        Bacteria          73 A          discrimina… dete… strong         NA
#> 2 Ala        Bacteria           3 G          acceptor_s… dete… strong         70
#> 3 Ala        Bacteria          70 U          acceptor_s… dete… strong          3
#> 4 Ala        Bacteria           4 G          acceptor_s… dete… weak           69
#> 5 Ala        Bacteria          69 C          acceptor_s… dete… weak            4
#> 6 Ala        Bacteria          20 G          d_arm       dete… weak           NA
#> # ℹ 4 more variables: pair_type <chr>, against_aars <chr>, universal <lgl>,
#> #   description <chr>
```
