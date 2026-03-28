# Summarize per-position modification frequency from modkit extract

Read a `mod_calls.tsv.gz` file from `modkit extract` and compute
per-position modification frequency (proportion of reads with a
non-canonical call code).

## Usage

``` r
summarize_mod_calls(path, refs = NULL, min_reads = 1)
```

## Arguments

- path:

  Path to a `mod_calls.tsv.gz` file.

- refs:

  Optional character vector of reference names to include.

- min_reads:

  Minimum number of reads at a position to retain. Default `1`.

## Value

A tibble with columns: `ref`, `pos`, `n_reads`, `n_mod`, `mod_freq`.

## Examples

``` r
path <- clover_example("ecoli/mod_calls.tsv.gz")
summarize_mod_calls(path)
#> # A tibble: 3 × 5
#>   ref                     pos n_reads n_mod mod_freq
#>   <chr>                 <dbl>   <int> <int>    <dbl>
#> 1 host-tRNA-Glu-TTC-1-1    34      15     3      0.2
#> 2 host-tRNA-Glu-TTC-1-1    46      15     3      0.2
#> 3 host-tRNA-Glu-TTC-1-1    55      15     0      0  
```
