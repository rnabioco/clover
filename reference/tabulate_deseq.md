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

tRNA-11

0.29

1.47 × 10⁻³

3.58 × 10⁻²

FALSE

tRNA-4

0.13

5.82 × 10⁻³

1.13 × 10⁻¹

TRUE

tRNA-18

−0.21

7.46 × 10⁻³

1.42 × 10⁻¹

FALSE

tRNA-3

0.44

1.66 × 10⁻²

1.52 × 10⁻¹

TRUE

tRNA-16

−1.43

3.60 × 10⁻²

1.66 × 10⁻¹

FALSE

tRNA-13

−0.60

4.53 × 10⁻²

1.43 × 10⁻¹

FALSE

tRNA-5

−0.83

4.59 × 10⁻²

1.56 × 10⁻¹

TRUE

tRNA-10

−0.83

4.97 × 10⁻²

9.13 × 10⁻²

TRUE

tRNA-17

−0.01

5.03 × 10⁻²

1.29 × 10⁻¹

FALSE

tRNA-9

1.46

5.24 × 10⁻²

1.69 × 10⁻¹

TRUE
