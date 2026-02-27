# Read modifications file

Read a TSV file with 4 columns: `ref`, `pos`, `mod_full`, `mod1`.

## Usage

``` r
read_mod_annotations(mods)
```

## Arguments

- mods:

  Path to a modifications TSV file.

## Value

A tibble.

## Examples

``` r
path <- clover_example("ecoli/mod_annotations.tsv")
read_mod_annotations(path)
#> # A tibble: 3 × 4
#>   ref                     pos mod_full          mod1 
#>   <chr>                 <dbl> <chr>             <chr>
#> 1 host-tRNA-Glu-TTC-1-1    34 5-methyluridine   m5U  
#> 2 host-tRNA-Glu-TTC-1-1    46 7-methylguanosine m7G  
#> 3 host-tRNA-Phe-GAA-1-1    34 5-methyluridine   m5U  
```
