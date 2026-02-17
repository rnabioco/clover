# Build a charging count matrix for differential charging analysis.

Create a count matrix with two columns per sample (charged and
uncharged), suitable for a DESeq2 interaction design that tests for
differences in charging ratios between conditions.

## Usage

``` r
charging_count_matrix(charging_data, min_count = 10)
```

## Arguments

- charging_data:

  A tibble of combined charging data with columns `tRNA`,
  `counts_charged`, `counts_uncharged`, and `sample_id` (as returned by
  [`read_charging_multi()`](https://rnabioco.github.io/clover/reference/read_charging_multi.md)).

- min_count:

  Minimum total count across all columns for a tRNA to be retained.
  Default `10`.

## Value

An integer matrix with tRNA names as row names. Column names are
`{sample_id}_charged` and `{sample_id}_uncharged`.

## Examples

``` r
if (FALSE) { # \dontrun{
charging <- read_charging_multi(paths)
mat <- charging_count_matrix(charging)
} # }
```
