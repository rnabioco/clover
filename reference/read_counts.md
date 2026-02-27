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
read_counts(clover_example("ecoli/counts.tsv"))
#> # A tibble: 3 × 2
#>   ref                   count
#>   <chr>                 <dbl>
#> 1 host-tRNA-Glu-TTC-1-1  1542
#> 2 host-tRNA-Ala-GGC-1-1    89
#> 3 host-tRNA-Phe-GAA-1-1   723
```
