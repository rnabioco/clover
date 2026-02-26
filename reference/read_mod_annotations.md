# Read modifications file

Read a TSV file with 4 columns: `ref`, `pos`, `mod_full`, `mod1`.

## Usage

``` r
read_mod_annotations(mods)
```

## Arguments

- mods:

  Path to a modifications TSV file.

## Value

A tibble.

## Examples

``` r
if (FALSE) { # \dontrun{
read_mod_annotations("modifications.tsv")
} # }
```
