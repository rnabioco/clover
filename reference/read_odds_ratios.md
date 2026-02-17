# Read an odds ratios file.

Read a per-position-pair odds ratio file produced by the tRNA sequencing
pipeline. These files contain pairwise modification co-occurrence
statistics.

## Usage

``` r
read_odds_ratios(path)
```

## Arguments

- path:

  Path to a `{sample}.odds_ratios.tsv.gz` file.

## Value

A tibble with columns including `pos1`, `pos2`, `odds_ratio`,
`log_odds_ratio`, `p_value`, and `total_obs`.

## Examples

``` r
if (FALSE) { # \dontrun{
or_data <- read_odds_ratios("sample1.odds_ratios.tsv.gz")
or_data
} # }
```
