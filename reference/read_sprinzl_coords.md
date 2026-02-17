# Read Sprinzl coordinates from a global coordinates TSV file.

Read a tRNAs-in-space global coordinates file and return a tibble with
standardized column names. These files map each position of each tRNA to
its Sprinzl numbering, structural region, and global alignment index.

## Usage

``` r
read_sprinzl_coords(path)
```

## Arguments

- path:

  Path to a global coordinates TSV file (may be gzipped).

## Value

A tibble with columns:

- `trna_id`: tRNA identifier (e.g., `nuc-tRNA-Ala-AGC-1-1`)

- `pos`: 1-based position in tRNA body

- `sprinzl_label`: Sprinzl position (character: "1", "20a", "47:e1",
  etc.)

- `global_index`: cross-tRNA universal alignment position

- `region`: structural region (acceptor-stem, D-stem, D-loop, etc.)

- `residue`: reference nucleotide

## Examples

``` r
path <- clover_example("sprinzl/sacCer_global_coords.tsv.gz")
coords <- read_sprinzl_coords(path)
coords
#> # A tibble: 19,967 × 6
#>    trna_id                pos sprinzl_label global_index region        residue
#>    <chr>                <dbl> <chr>                <dbl> <chr>         <chr>  
#>  1 nuc-tRNA-Ala-AGC-1-1     1 1                        1 acceptor-stem G      
#>  2 nuc-tRNA-Ala-AGC-1-1     2 2                        2 acceptor-stem G      
#>  3 nuc-tRNA-Ala-AGC-1-1     3 3                        3 acceptor-stem G      
#>  4 nuc-tRNA-Ala-AGC-1-1     4 4                        4 acceptor-stem C      
#>  5 nuc-tRNA-Ala-AGC-1-1     5 5                        5 acceptor-stem G      
#>  6 nuc-tRNA-Ala-AGC-1-1     6 6                        6 acceptor-stem T      
#>  7 nuc-tRNA-Ala-AGC-1-1     7 7                        8 acceptor-stem G      
#>  8 nuc-tRNA-Ala-AGC-1-1     8 8                        9 unknown       T      
#>  9 nuc-tRNA-Ala-AGC-1-1     9 9                       10 unknown       G      
#> 10 nuc-tRNA-Ala-AGC-1-1    10 10                      11 D-stem        G      
#> # ℹ 19,957 more rows
```
