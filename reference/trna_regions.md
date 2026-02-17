# tRNA structural regions.

Return a named list mapping canonical tRNA structural region names to
integer vectors of Sprinzl position numbers.

## Usage

``` r
trna_regions()
```

## Value

A named list of integer vectors.

## Examples

``` r
trna_regions()
#> $acceptor_stem
#>  [1]  1  2  3  4  5  6  7 66 67 68 69 70 71 72
#> 
#> $d_arm
#>  [1]  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25
#> 
#> $anticodon_stem
#>  [1] 26 27 28 29 30 31 39 40 41 42 43 44
#> 
#> $anticodon_loop
#> [1] 32 33 34 35 36 37 38
#> 
#> $variable_loop
#> [1] 44 45 46 47 48
#> 
#> $t_arm
#>  [1] 49 50 51 52 53 54 55 56 57 58 59 60 61 62 63 64 65
#> 
#> $discriminator
#> [1] 73
#> 
#> $cca
#> [1] 74 75 76
#> 
```
