# Tabulate top DESeq2 differential expression results.

Create a formatted gt table of the top significant tRNAs from
[`tidy_deseq_results()`](https://rnabioco.github.io/clover/reference/tidy_deseq_results.md)
output, sorted by p-value.

## Usage

``` r
tabulate_deseq(data, lab_col = "ref", n = 10)
```

## Arguments

- data:

  A tibble from
  [`tidy_deseq_results()`](https://rnabioco.github.io/clover/reference/tidy_deseq_results.md)
  with at least `log2FoldChange`, `pvalue`, `padj`, and `significant`
  columns.

- lab_col:

  Column name (string) used for row labels. Default `"ref"`.

- n:

  Maximum number of rows to display. Default `10`.

## Value

A `gt_tbl` object.

## Examples

``` r
res <- tibble::tibble(
  ref = paste0("tRNA-", 1:20),
  log2FoldChange = rnorm(20),
  pvalue = runif(20, 0, 0.1),
  padj = runif(20, 0, 0.2),
  significant = c(rep(TRUE, 10), rep(FALSE, 10))
)
if (requireNamespace("gt", quietly = TRUE)) {
  tabulate_deseq(res)
}


  

ref
```

log2 FC

p-value

Adjusted p-value

significant

tRNA-9

−0.83

1.47 × 10⁻³

3.58 × 10⁻²

TRUE

tRNA-2

0.44

5.82 × 10⁻³

1.13 × 10⁻¹

TRUE

tRNA-16

−0.01

7.46 × 10⁻³

1.42 × 10⁻¹

FALSE

tRNA-1

−0.08

1.66 × 10⁻²

1.52 × 10⁻¹

TRUE

tRNA-14

0.15

3.60 × 10⁻²

1.66 × 10⁻¹

FALSE

tRNA-11

−0.48

4.53 × 10⁻²

1.43 × 10⁻¹

FALSE

tRNA-3

0.13

4.59 × 10⁻²

1.56 × 10⁻¹

TRUE

tRNA-8

1.46

4.97 × 10⁻²

9.13 × 10⁻²

TRUE

tRNA-15

−1.43

5.03 × 10⁻²

1.29 × 10⁻¹

FALSE

tRNA-7

−0.75

5.24 × 10⁻²

1.69 × 10⁻¹

TRUE
