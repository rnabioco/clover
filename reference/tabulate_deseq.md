# Tabulate top DESeq2 differential expression results.

Create a formatted gt table of the top significant tRNAs from
[`tidy_deseq_results()`](https://rnabioco.github.io/clover/reference/tidy_deseq_results.md)
output, sorted by p-value.

## Usage

``` r
tabulate_deseq(data, lab_col = "tRNA", n = 10)
```

## Arguments

- data:

  A tibble from
  [`tidy_deseq_results()`](https://rnabioco.github.io/clover/reference/tidy_deseq_results.md)
  with at least `log2FoldChange`, `pvalue`, `padj`, and `significant`
  columns.

- lab_col:

  Column name (string) used for row labels. Default `"tRNA"`.

- n:

  Maximum number of rows to display. Default `10`.

## Value

A `gt_tbl` object.

## Examples

``` r
res <- tibble::tibble(
  tRNA = paste0("tRNA-", 1:20),
  log2FoldChange = rnorm(20),
  pvalue = runif(20, 0, 0.1),
  padj = runif(20, 0, 0.2),
  significant = c(rep(TRUE, 10), rep(FALSE, 10))
)
if (requireNamespace("gt", quietly = TRUE)) {
  tabulate_deseq(res)
}


  

tRNA
```

log2 FC

p-value

Adjusted p-value

significant

tRNA-5

1.61

1.48 × 10⁻³

1.61 × 10⁻¹

TRUE

tRNA-18

−1.14

2.87 × 10⁻³

9.18 × 10⁻²

FALSE

tRNA-15

−1.80

3.09 × 10⁻³

1.42 × 10⁻¹

FALSE

tRNA-2

0.12

1.28 × 10⁻²

7.57 × 10⁻²

TRUE

tRNA-17

−0.24

1.67 × 10⁻²

1.85 × 10⁻¹

FALSE

tRNA-13

1.20

1.70 × 10⁻²

6.55 × 10⁻²

FALSE

tRNA-3

−0.28

2.61 × 10⁻²

3.49 × 10⁻²

TRUE

tRNA-4

0.50

2.88 × 10⁻²

1.21 × 10⁻¹

TRUE

tRNA-12

−0.92

3.61 × 10⁻²

1.28 × 10⁻¹

FALSE

tRNA-7

0.39

4.14 × 10⁻²

1.47 × 10⁻¹

TRUE
