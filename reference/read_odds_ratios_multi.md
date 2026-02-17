# Read odds ratio files for multiple samples.

Read and combine odds ratio files from multiple samples into a single
tibble with a `sample_id` column.

## Usage

``` r
read_odds_ratios_multi(paths)
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
paths <- c(sample1 = "s1.odds_ratios.tsv.gz",
           sample2 = "s2.odds_ratios.tsv.gz")
or_data <- read_odds_ratios_multi(paths)
} # }
```
