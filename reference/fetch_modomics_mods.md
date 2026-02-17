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
if (FALSE) { # \dontrun{
fa <- clover_example("yeast/trna-ref.fa.gz")
mods <- fetch_modomics_mods(fa, "Saccharomyces cerevisiae")
mods
} # }
```
