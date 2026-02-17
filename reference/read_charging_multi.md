# Read charging CPM files for multiple samples.

Read and combine charging CPM files from multiple samples into a single
tibble with a `sample_id` column.

## Usage

``` r
read_charging_multi(paths)
```

## Arguments

- paths:

  A named character vector of file paths. Names are used as sample
  identifiers.

## Value

A tibble with all samples combined and a `sample_id` column.

## Examples

``` r
if (FALSE) { # \dontrun{
paths <- c(sample1 = "s1.charging.cpm.tsv.gz",
           sample2 = "s2.charging.cpm.tsv.gz")
charging <- read_charging_multi(paths)
} # }
```
