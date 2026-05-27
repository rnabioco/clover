# Retrieve canonical tRNA tertiary contacts

Returns the canonical cloverleaf tertiary interactions that stabilize
the L-shape and elbow of the tRNA fold: four base triples and five
tertiary base pairs. Numbering follows the Sprinzl convention. Modified
bases (m7G46, Psi55, m1A58) are encoded as their unmodified equivalents
(G, U, A) so callers can compare against in-vitro transcripts directly.

## Usage

``` r
tertiary_contacts(organism)
```

## Arguments

- organism:

  Character string specifying the organism name (e.g.,
  `"Escherichia coli"`, `"Saccharomyces cerevisiae"`). Use
  [`identity_organisms()`](https://rnabioco.github.io/clover/reference/identity_organisms.md)
  to list supported organisms.

## Value

A tibble with one row per contact and columns:

- `contact_id`: short label (e.g., `"U8-A14-A21"`, `"G15-C48"`)

- `contact_type`: `"triple"` or `"pair"`

- `domain`: `"Bacteria"` or `"Eukarya"`

- `pos1`, `pos2`, `pos3`: Sprinzl positions (`pos3 = NA` for pairs)

- `nuc1`, `nuc2`, `nuc3`: canonical bases (`nuc3 = NA` for pairs)

- `interaction`: geometry label (e.g., `"reverse-Hoogsteen"`,
  `"trans-WC"`, `"Levitt"`)

- `universal`: `TRUE` for contacts conserved across all domains

- `description`: human-readable description

## Details

Sources: Westhof & Auffinger (2012) and Giege & Eriani (2023). These
tertiary contacts are not aaRS identity elements per se but constrain
the tRNA fold and may indirectly affect aminoacylation specificity. Use
[`identity_elements()`](https://rnabioco.github.io/clover/reference/identity_elements.md)
for aaRS recognition determinants.

## Limitations

Selenocysteine tRNA (SeC, anticodon UCA) is not currently supported. Its
non-canonical 90-nt fold with an extended variable arm does not share
the canonical tertiary contact set bundled here, and SeC Sprinzl
coordinates are not available.

## References

Westhof E, Auffinger P (2012). "tRNA structure." *Encyclopedia of Life
Sciences*. John Wiley & Sons.

Giege R, Eriani G (2023). "The tRNA identity landscape for
aminoacylation and beyond." *Nucleic Acids Research*, 51(4), 1528–1570.
[doi:10.1093/nar/gkad007](https://doi.org/10.1093/nar/gkad007)

## Examples

``` r
tertiary_contacts("Escherichia coli")
#> # A tibble: 9 × 12
#>   contact_id contact_type domain  pos1  pos2  pos3 nuc1  nuc2  nuc3  interaction
#>   <chr>      <chr>        <chr>  <int> <int> <int> <chr> <chr> <chr> <chr>      
#> 1 U8-A14-A21 triple       Bacte…     8    14    21 U     A     A     reverse-Ho…
#> 2 A9-U12-A23 triple       Bacte…     9    12    23 A     U     A     trans-WC   
#> 3 G10-C25-G… triple       Bacte…    10    25    45 G     C     G     trans-WC   
#> 4 C13-G22-G… triple       Bacte…    13    22    46 C     G     G     trans-WC   
#> 5 G15-C48    pair         Bacte…    15    48    NA G     C     NA    Levitt     
#> 6 G18-U55    pair         Bacte…    18    55    NA G     U     NA    trans-WC   
#> 7 G19-C56    pair         Bacte…    19    56    NA G     C     NA    WC         
#> 8 A26-G44    pair         Bacte…    26    44    NA A     G     NA    non-canoni…
#> 9 U54-A58    pair         Bacte…    54    58    NA U     A     NA    reverse-Ho…
#> # ℹ 2 more variables: universal <lgl>, description <chr>

# Just the four classical triples
tc <- tertiary_contacts("Escherichia coli")
tc[tc$contact_type == "triple", ]
#> # A tibble: 4 × 12
#>   contact_id contact_type domain  pos1  pos2  pos3 nuc1  nuc2  nuc3  interaction
#>   <chr>      <chr>        <chr>  <int> <int> <int> <chr> <chr> <chr> <chr>      
#> 1 U8-A14-A21 triple       Bacte…     8    14    21 U     A     A     reverse-Ho…
#> 2 A9-U12-A23 triple       Bacte…     9    12    23 A     U     A     trans-WC   
#> 3 G10-C25-G… triple       Bacte…    10    25    45 G     C     G     trans-WC   
#> 4 C13-G22-G… triple       Bacte…    13    22    46 C     G     G     trans-WC   
#> # ℹ 2 more variables: universal <lgl>, description <chr>
```
