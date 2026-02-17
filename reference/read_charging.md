# Read a charging CPM file.

Read a per-tRNA charging CPM file produced by the tRNA sequencing
pipeline. These files contain counts and CPM values for charged and
uncharged reads.

## Usage

``` r
read_charging(path)
```

## Arguments

- path:

  Path to a `{sample}.charging.cpm.tsv.gz` file.

## Value

A tibble with columns including `tRNA`, `counts_charged`,
`counts_uncharged`, `cpm_charged`, `cpm_uncharged`, and `total_count`.

## Examples

``` r
if (FALSE) { # \dontrun{
charging <- read_charging("sample1.charging.cpm.tsv.gz")
charging
} # }
```
