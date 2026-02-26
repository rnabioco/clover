# Read FASTA reference

Read FASTA reference

## Usage

``` r
read_fasta(fa)
```

## Arguments

- fa:

  Path to a FASTA file.

## Value

A
[Biostrings::DNAStringSet](https://rdrr.io/pkg/Biostrings/man/XStringSet-class.html).

## Examples

``` r
fa <- clover_example("ecoli/trna_only.fa.gz")
read_fasta(fa)
#> DNAStringSet object of length 95:
#>      width seq                                              names               
#>  [1]    77 AGGCTTGTAGCTCAGGTGGTTAG...TCAAGTCCACTCAGGCCTACCA host-tRNA-Ile-GAT...
#>  [2]    76 GGGGCTATAGCTCAGCTGGGAGA...TCGATCCCGCATAGCTCCACCA host-tRNA-Ala-TGC...
#>  [3]    77 GGAGCGGTAGTTCAGTCGGTTAG...TCGAGTCCCGTCCGTTCCGCCA host-tRNA-Asp-GTC...
#>  [4]    77 GGAGCGGTAGTTCAGTCGGTTAG...TCGAGTCCCGTCCGTTCCGCCA host-tRNA-Asp-GTC...
#>  [5]    76 GCCGATATAGCTCAGTTGGTAGA...TCGACTCCTATTATCGGCACCA host-tRNA-Thr-CGT...
#>  ...   ... ...
#> [91]    78 CTCCGTGTAGCTCAGTTTGGTAG...CAAATCCTTGTATGGAGAGCCA phage-tRNA-Pro-TGG
#> [92]    90 GGAGGCGTGGCAGAGTGGTTTAA...TCAAATCCTATCGCCTCCGCCA phage-tRNA-Ser-TGA
#> [93]    76 GCTGATTTAGCTCAGTAGGTAGA...TCGATTCCGTCAATCAGCACCA phage-tRNA-Thr-TGT
#> [94]    77 GGCCCTGTAGCTCAATGGTTAGC...TCAAATCTGGTCTGGGTCACCA phage-tRNA-Ile2-CAT
#> [95]    75 GTCCCGCTGGTGTAATGGATAGC...TCGATCCCAGGGCGGGATACCA phage-tRNA-Arg-TCT
```
