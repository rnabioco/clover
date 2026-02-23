# Tidy DESeq2 results into a tibble.

Extract results from a DESeq2 analysis and return a tidy tibble with
tRNA identifiers and significance flags.

## Usage

``` r
tidy_deseq_results(dds, contrast, padj_cutoff = 0.05)
```

## Arguments

- dds:

  A `DESeqDataSet` object (from
  [`run_deseq()`](https://rnabioco.github.io/clover/reference/run_deseq.md)).

- contrast:

  A contrast specification: either a character vector of length 3 (e.g.,
  `c("condition", "mutant", "wildtype")`) or a list for
  coefficient-based contrasts.

- padj_cutoff:

  Adjusted p-value threshold for significance. Default `0.05`.

## Value

A tibble with columns: `ref`, `log2FoldChange`, `lfcSE`, `pvalue`,
`padj`, and `significant` (logical).

## Examples

``` r
# \donttest{
se <- create_clover(clover_example("ecoli/config.yaml"))
counts <- SummarizedExperiment::assay(se, "counts")
coldata <- as.data.frame(SummarizedExperiment::colData(se))
coldata$condition <- ifelse(
  grepl("ctl", coldata$sample_id), "ctl", "inf"
)
dds <- run_deseq(counts, coldata, design = ~condition)
#> Warning: some variables in design formula are characters, converting to factors
#> estimating size factors
#> estimating dispersions
#> gene-wise dispersion estimates
#> mean-dispersion relationship
#> final dispersion estimates
#> fitting model and testing
tidy_deseq_results(dds, contrast = c("condition", "inf", "ctl"))
#> # A tibble: 190 × 6
#>    ref                             log2FoldChange lfcSE pvalue  padj significant
#>    <chr>                                    <dbl> <dbl>  <dbl> <dbl> <lgl>      
#>  1 host-tRNA-Ala-GGC-1-1                  -0.236  0.393  0.548 0.822 FALSE      
#>  2 host-tRNA-Ala-GGC-1-1-uncharged         0.163  0.352  0.643 0.857 FALSE      
#>  3 host-tRNA-Ala-GGC-1-2                   0.145  0.434  0.739 0.896 FALSE      
#>  4 host-tRNA-Ala-GGC-1-2-uncharged         0.166  0.360  0.645 0.857 FALSE      
#>  5 host-tRNA-Ala-TGC-1-1                  -0.424  0.350  0.225 0.769 FALSE      
#>  6 host-tRNA-Ala-TGC-1-1-uncharged         0.0849 0.267  0.750 0.900 FALSE      
#>  7 host-tRNA-Ala-TGC-1-2                  -0.314  0.395  0.427 0.790 FALSE      
#>  8 host-tRNA-Ala-TGC-1-2-uncharged         0.260  0.472  0.581 0.832 FALSE      
#>  9 host-tRNA-Ala-TGC-1-3                   0.0823 0.427  0.847 0.922 FALSE      
#> 10 host-tRNA-Ala-TGC-1-3-uncharged         0.112  0.294  0.703 0.896 FALSE      
#> # ℹ 180 more rows
# }
```
