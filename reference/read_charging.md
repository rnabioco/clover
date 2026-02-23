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

A tibble with columns including `ref`, `counts_charged`,
`counts_uncharged`, `cpm_charged`, `cpm_uncharged`, and `total_count`.

## Examples

``` r
path <- clover_example(
  "ecoli/summary/tables/wt-15-ctl-01/wt-15-ctl-01.charging.cpm.tsv.gz"
)
read_charging(path)
#> # A tibble: 190 × 5
#>    ref                 counts_charged counts_uncharged cpm_charged cpm_uncharged
#>    <chr>                        <dbl>            <dbl>       <dbl>         <dbl>
#>  1 host-tRNA-Ala-GGC-…            174              124        270.          193.
#>  2 host-tRNA-Ala-GGC-…           1173             3220       1824.         5006.
#>  3 host-tRNA-Ala-GGC-…            818              807       1272.         1255.
#>  4 host-tRNA-Ala-GGC-…           1205             3186       1873.         4953.
#>  5 host-tRNA-Ala-TGC-…            351              167        546.          260.
#>  6 host-tRNA-Ala-TGC-…           1390             2395       2161.         3723.
#>  7 host-tRNA-Ala-TGC-…            615              214        956.          333.
#>  8 host-tRNA-Ala-TGC-…           4739            13203       7367.        20525.
#>  9 host-tRNA-Ala-TGC-…           1167             1162       1814.         1806.
#> 10 host-tRNA-Ala-TGC-…           1321             2330       2054.         3622.
#> # ℹ 180 more rows
```
