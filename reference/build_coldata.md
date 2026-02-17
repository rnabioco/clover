# Build column metadata for DESeq2.

Create a `data.frame` of sample metadata suitable for
[`DESeq2::DESeqDataSetFromMatrix()`](https://rdrr.io/pkg/DESeq2/man/DESeqDataSet.html).
For charging count matrices, the `charge_status` factor is added
automatically.

## Usage

``` r
build_coldata(count_matrix, sample_info = NULL)
```

## Arguments

- count_matrix:

  A count matrix from
  [`abundance_count_matrix()`](https://rnabioco.github.io/clover/reference/abundance_count_matrix.md)
  or
  [`charging_count_matrix()`](https://rnabioco.github.io/clover/reference/charging_count_matrix.md).

- sample_info:

  An optional data frame with a `sample_id` column and additional
  experimental factor columns (e.g., `condition`, `replicate`). If
  `NULL`, a minimal data frame is created from column names.

## Value

A data frame with row names matching `colnames(count_matrix)`.

## Examples

``` r
if (FALSE) { # \dontrun{
mat <- abundance_count_matrix(charging)
sample_info <- data.frame(
  sample_id = c("wt_1", "wt_2", "mut_1", "mut_2"),
  condition = c("wt", "wt", "mut", "mut")
)
coldata <- build_coldata(mat, sample_info)
} # }
```
