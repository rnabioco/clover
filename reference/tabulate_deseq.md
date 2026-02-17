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

tRNA-7

−0.40

1.95 × 10⁻³

1.46 × 10⁻²

TRUE

tRNA-8

−0.82

8.72 × 10⁻³

1.77 × 10⁻¹

TRUE

tRNA-17

1.36

1.02 × 10⁻²

1.10 × 10⁻¹

FALSE

tRNA-12

0.38

1.15 × 10⁻²

1.05 × 10⁻¹

FALSE

tRNA-11

0.13

1.24 × 10⁻²

7.35 × 10⁻²

FALSE

tRNA-3

2.56

1.36 × 10⁻²

1.42 × 10⁻¹

TRUE

tRNA-5

1.14

2.63 × 10⁻²

1.32 × 10⁻¹

TRUE

tRNA-13

1.14

3.67 × 10⁻²

9.40 × 10⁻²

FALSE

tRNA-4

1.06

3.76 × 10⁻²

1.94 × 10⁻¹

TRUE

tRNA-14

1.24

4.71 × 10⁻²

1.97 × 10⁻²

FALSE
