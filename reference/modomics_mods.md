# Map cached MODOMICS modifications onto reference sequences

Uses bundled MODOMICS data to map known tRNA modification positions onto
user-provided reference sequences using pairwise alignment. No internet
connection is required for organisms included in the package (see
[`modomics_organisms()`](https://rnabioco.github.io/clover/reference/modomics_organisms.md)).
For other organisms, falls back to
[`fetch_modomics_mods()`](https://rnabioco.github.io/clover/reference/fetch_modomics_mods.md).

## Usage

``` r
modomics_mods(fasta, organism, min_identity = 0.7)
```

## Arguments

- fasta:

  Path to a FASTA file or a
  [Biostrings::DNAStringSet](https://rdrr.io/pkg/Biostrings/man/XStringSet-class.html)
  object containing reference tRNA sequences.

- organism:

  Character string specifying the organism name as used in MODOMICS
  (e.g., `"Saccharomyces cerevisiae"`, `"Escherichia coli"`).

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
fa <- clover_example("ecoli/trna_only.fa.gz")
modomics_mods(fa, "Escherichia coli")
#> Processing 182 MODOMICS sequences.
#> Matching MODOMICS sequences to reference FASTA.
#> Found 701 modification annotations.
#> # A tibble: 701 × 4
#>    ref                     pos mod_full                 mod1 
#>    <chr>                 <int> <chr>                    <chr>
#>  1 host-tRNA-Ala-TGC-1-1    17 dihydrouridine           D    
#>  2 host-tRNA-Ala-TGC-1-1    34 uridine 5-oxyacetic acid cmo5U
#>  3 host-tRNA-Ala-TGC-1-1    46 7-methylguanosine        m7G  
#>  4 host-tRNA-Ala-TGC-1-1    54 5-methyluridine          m5U  
#>  5 host-tRNA-Ala-TGC-1-1    55 pseudouridine            Y    
#>  6 host-tRNA-Ala-GGC-1-1    17 dihydrouridine           D    
#>  7 host-tRNA-Ala-GGC-1-1    46 7-methylguanosine        m7G  
#>  8 host-tRNA-Ala-GGC-1-1    54 5-methyluridine          m5U  
#>  9 host-tRNA-Ala-GGC-1-1    55 pseudouridine            Y    
#> 10 host-tRNA-Arg-ACG-1-1     8 4-thiouridine            s4U  
#> # ℹ 691 more rows
# }
```
