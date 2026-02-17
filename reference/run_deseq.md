# Run DESeq2 differential analysis.

Wrapper around
[`DESeq2::DESeqDataSetFromMatrix()`](https://rdrr.io/pkg/DESeq2/man/DESeqDataSet.html)
and [`DESeq2::DESeq()`](https://rdrr.io/pkg/DESeq2/man/DESeq.html) for
tRNA abundance or charging data.

## Usage

``` r
run_deseq(count_matrix, coldata, design, ...)
```

## Arguments

- count_matrix:

  An integer count matrix (from
  [`abundance_count_matrix()`](https://rnabioco.github.io/clover/reference/abundance_count_matrix.md)
  or
  [`charging_count_matrix()`](https://rnabioco.github.io/clover/reference/charging_count_matrix.md)).

- coldata:

  A data frame of sample metadata (from
  [`build_coldata()`](https://rnabioco.github.io/clover/reference/build_coldata.md)).

- design:

  A formula specifying the design (e.g., `~ condition` or
  `~ condition + charge_status + condition:charge_status`).

- ...:

  Additional arguments passed to
  [`DESeq2::DESeq()`](https://rdrr.io/pkg/DESeq2/man/DESeq.html).

## Value

A `DESeqDataSet` object.

## Examples

``` r
if (FALSE) { # \dontrun{
dds <- run_deseq(mat, coldata, design = ~ condition)
} # }
```
