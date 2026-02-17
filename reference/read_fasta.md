# Read FASTA reference

Read FASTA reference

## Usage

``` r
read_fasta(fa)
```

## Arguments

- fa:

  path to fasta file

## Value

A
[Biostrings::DNAStringSet](https://rdrr.io/pkg/Biostrings/man/XStringSet-class.html).

## Examples

``` r
fa <- clover_example("ecoli/validated.fa.gz")
read_fasta(fa)
#> DNAStringSet object of length 190:
#>       width seq                                             names               
#>   [1]   141 CCTAAGAGCAAGAAGAAGCCTG...CAACCTTGCCTTAAAAAAAAAA host-tRNA-Ile-GAT...
#>   [2]   141 CCTAAGAGCAAGAAGAAGCCTG...TGGAAGGTAGGCAAAAAAAAAA host-tRNA-Ile-GAT...
#>   [3]   140 CCTAAGAGCAAGAAGAAGCCTG...CAACCTTGCCTTAAAAAAAAAA host-tRNA-Ala-TGC...
#>   [4]   140 CCTAAGAGCAAGAAGAAGCCTG...TGGAAGGTAGGCAAAAAAAAAA host-tRNA-Ala-TGC...
#>   [5]   141 CCTAAGAGCAAGAAGAAGCCTG...CAACCTTGCCTTAAAAAAAAAA host-tRNA-Asp-GTC...
#>   ...   ... ...
#> [186]   140 CCTAAGAGCAAGAAGAAGCCTG...TGGAAGGTAGGCAAAAAAAAAA phage-tRNA-Thr-TG...
#> [187]   141 CCTAAGAGCAAGAAGAAGCCTG...CAACCTTGCCTTAAAAAAAAAA phage-tRNA-Ile2-CAT
#> [188]   141 CCTAAGAGCAAGAAGAAGCCTG...TGGAAGGTAGGCAAAAAAAAAA phage-tRNA-Ile2-C...
#> [189]   139 CCTAAGAGCAAGAAGAAGCCTG...CAACCTTGCCTTAAAAAAAAAA phage-tRNA-Arg-TCT
#> [190]   139 CCTAAGAGCAAGAAGAAGCCTG...TGGAAGGTAGGCAAAAAAAAAA phage-tRNA-Arg-TC...
```
