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
if (FALSE) { # \dontrun{
fa <- clover_example("ecoli/validated.fa.gz")
mods <- modomics_mods(fa, "Escherichia coli")
mods
} # }
```
