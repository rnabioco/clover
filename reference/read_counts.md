# Read a counts TSV file.

Read a per-tRNA counts file produced by the tRNA sequencing pipeline.

## Usage

``` r
read_counts(path)
```

## Arguments

- path:

  Path to a counts TSV file (may be gzipped).

## Value

A tibble.

## Examples

``` r
if (FALSE) { # \dontrun{
counts <- read_counts("sample1.counts.tsv.gz")
} # }
```
