# Compute pairwise modification co-occurrence odds ratios.

Given a modkit `mod_calls.tsv.gz` file, build a per-read binary
modification matrix for each tRNA and compute Fisher's exact test for
each pair of positions. This is computationally expensive for full
datasets; use the `refs` parameter to restrict to specific tRNAs.

## Usage

``` r
compute_odds_ratios(mod_calls_path, refs = NULL, min_reads = 10)
```

## Arguments

- mod_calls_path:

  Path to a `{sample}.mod_calls.tsv.gz` file from modkit.

- refs:

  Optional character vector of reference names to include. If `NULL`,
  all references are processed.

- min_reads:

  Minimum number of reads required for a tRNA to be included. Default
  `10`.

## Value

A tibble with columns: `ref`, `pos1`, `pos2`, `odds_ratio`,
`log_odds_ratio`, `p_value`, `total_obs`.
