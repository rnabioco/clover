# List available tRNA structures for an organism

Returns the names of tRNAs for which cloverleaf structure SVGs are
bundled with the package for the given organism.

## Usage

``` r
structure_trnas(organism)
```

## Arguments

- organism:

  Character string specifying the organism name (e.g.,
  `"Escherichia coli"`). Use
  [`structure_organisms()`](https://rnabioco.github.io/clover/reference/structure_organisms.md)
  to list available organisms.

## Value

A character vector of tRNA names.

## Examples

``` r
if (FALSE) { # \dontrun{
structure_trnas("Escherichia coli")
} # }
```
