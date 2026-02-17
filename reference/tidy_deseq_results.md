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
if (FALSE) { # \dontrun{
res <- tidy_deseq_results(dds, contrast = c("condition", "mut", "wt"))
res
} # }
```
