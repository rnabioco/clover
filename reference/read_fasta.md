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
fa <- clover_example("yeast/trna-ref.fa.gz")
read_fasta(fa)
#> DNAStringSet object of length 297:
#>       width seq                                             names               
#>   [1]   130 CCTAAGAGCAAGAAGAAGCCTG...CTTGCTCTTAGGAAAAAAAAAA nuc-tRNA-Ala-AGC-1-1
#>   [2]   130 CCTAAGAGCAAGAAGAAGCCTG...CTTGCTCTTAGGAAAAAAAAAA nuc-tRNA-Ala-AGC-...
#>   [3]   130 CCTAAGAGCAAGAAGAAGCCTG...CTTGCTCTTAGGAAAAAAAAAA nuc-tRNA-Ala-AGC-...
#>   [4]   130 CCTAAGAGCAAGAAGAAGCCTG...CTTGCTCTTAGGAAAAAAAAAA nuc-tRNA-Ala-AGC-1-2
#>   [5]   130 CCTAAGAGCAAGAAGAAGCCTG...CTTGCTCTTAGGAAAAAAAAAA nuc-tRNA-Ala-AGC-1-3
#>   ...   ... ...
#> [293]   131 CCTAAGAGCAAGAAGAAGCCTG...CTTGCTCTTAGGAAAAAAAAAA mito-tRNA-Met-CAU
#> [294]   129 CCTAAGAGCAAGAAGAAGCCTG...CTTGCTCTTAGGAAAAAAAAAA mito-tRNA-Phe-GAA
#> [295]   129 CCTAAGAGCAAGAAGAAGCCTG...CTTGCTCTTAGGAAAAAAAAAA mito-tRNA-Thr-UAG
#> [296]   130 CCTAAGAGCAAGAAGAAGCCTG...CTTGCTCTTAGGAAAAAAAAAA mito-tRNA-Val-UAC
#> [297]   130 CCTAAGAGCAAGAAGAAGCCTG...CTTGCTCTTAGGAAAAAAAAAA mito-tRNA-fMet-CAU
```
