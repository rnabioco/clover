# Fetch tRNA modification annotations from MODOMICS

Downloads tRNA modification data from the
[MODOMICS](https://genesilico.pl/modomics/) database and maps
modification positions onto user-provided reference sequences using
pairwise alignment.

## Usage

``` r
fetch_modomics_mods(fasta, organism, cache_dir = NULL, min_identity = 0.7)
```

## Arguments

- fasta:

  Path to a FASTA file or a
  [Biostrings::DNAStringSet](https://rdrr.io/pkg/Biostrings/man/XStringSet-class.html)
  object containing reference tRNA sequences.

- organism:

  Character string specifying the organism name as used in MODOMICS
  (e.g., `"Saccharomyces cerevisiae"`, `"Escherichia coli"`).

- cache_dir:

  Optional directory path for caching API responses as RDS files. If
  `NULL` (default), no caching is performed.

- min_identity:

  Minimum alignment identity (0–1) required to accept a match between a
  MODOMICS sequence and a reference sequence. Default `0.7`.

## Value

A tibble with columns:

- `ref`: reference sequence name from the FASTA

- `pos`: 1-based position in the reference sequence

- `mod_full`: full modification name (e.g., "1-methyladenosine")

- `mod1`: short modification name (e.g., "m1A")

## Examples

``` r
# \donttest{
fa <- clover_example("ecoli/validated.fa.gz")
mods <- fetch_modomics_mods(fa, "Escherichia coli")
#> Fetching MODOMICS modification dictionary.
#> Fetching MODOMICS tRNA sequences for "Escherichia coli".
#> Processing 182 MODOMICS sequences.
#> Matching MODOMICS sequences to reference FASTA.
#> Found 695 modification annotations.
mods
#> # A tibble: 695 × 4
#>    ref                     pos mod_full                 mod1 
#>    <chr>                 <int> <chr>                    <chr>
#>  1 host-tRNA-Ala-TGC-1-1    41 dihydrouridine           D    
#>  2 host-tRNA-Ala-TGC-1-1    58 uridine 5-oxyacetic acid cmo5U
#>  3 host-tRNA-Ala-TGC-1-1    70 7-methylguanosine        m7G  
#>  4 host-tRNA-Ala-TGC-1-1    78 5-methyluridine          m5U  
#>  5 host-tRNA-Ala-TGC-1-1    79 pseudouridine            Y    
#>  6 host-tRNA-Ala-GGC-1-1    41 dihydrouridine           D    
#>  7 host-tRNA-Ala-GGC-1-1    70 7-methylguanosine        m7G  
#>  8 host-tRNA-Ala-GGC-1-1    78 5-methyluridine          m5U  
#>  9 host-tRNA-Ala-GGC-1-1    79 pseudouridine            Y    
#> 10 host-tRNA-Arg-ACG-1-1    32 4-thiouridine            s4U  
#> # ℹ 685 more rows
# }
```
