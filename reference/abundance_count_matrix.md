# Build a tRNA abundance count matrix from charging data.

Pivot charging data into an integer matrix suitable for DESeq2, where
each row is a tRNA and each column is a sample. Values are total
abundance (charged + uncharged counts).

## Usage

``` r
abundance_count_matrix(charging_data, min_count = 10)
```

## Arguments

- charging_data:

  A tibble of combined charging data with columns `ref`,
  `counts_charged`, `counts_uncharged`, and `sample_id` (as returned by
  [`read_charging_multi()`](https://rnabioco.github.io/clover/reference/read_charging_multi.md)).

- min_count:

  Minimum total count across all samples for a tRNA to be retained.
  Default `10`.

## Value

An integer matrix with tRNA names as row names and sample IDs as column
names.

## Examples

``` r
if (FALSE) { # \dontrun{
charging <- read_charging_multi(paths)
mat <- abundance_count_matrix(charging)
} # }
```
