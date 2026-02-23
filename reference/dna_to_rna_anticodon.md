# Convert DNA anticodon to RNA in tRNA name strings

Replaces T with U in the anticodon portion of tRNA identifiers, e.g.,
`"tRNA-Glu-TTC-1-1"` becomes `"tRNA-Glu-UUC-1-1"`.

## Usage

``` r
dna_to_rna_anticodon(x)
```

## Arguments

- x:

  Character vector of tRNA name strings.

## Value

Character vector with DNA anticodons converted to RNA.

## Examples

``` r
dna_to_rna_anticodon("tRNA-Glu-TTC-1-1")
#> [1] "tRNA-Glu-UUC-1-1"
```
