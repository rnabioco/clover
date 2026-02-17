test_that("tabulate_deseq returns a gt_tbl", {
  skip_if_not_installed("gt")
  res <- tibble::tibble(
    ref = paste0("tRNA-", 1:20),
    log2FoldChange = rnorm(20),
    pvalue = runif(20, 0, 0.1),
    padj = runif(20, 0, 0.2),
    significant = c(rep(TRUE, 10), rep(FALSE, 10))
  )
  tbl <- tabulate_deseq(res)
  expect_s3_class(tbl, "gt_tbl")
})

test_that("tabulate_deseq respects n parameter", {
  skip_if_not_installed("gt")
  res <- tibble::tibble(
    ref = paste0("tRNA-", 1:20),
    log2FoldChange = rnorm(20),
    pvalue = runif(20, 0, 0.1),
    padj = runif(20, 0, 0.2),
    significant = rep(TRUE, 20)
  )
  tbl <- tabulate_deseq(res, n = 5)
  expect_equal(nrow(tbl[["_data"]]), 5)
})

test_that("tabulate_deseq filters NA padj rows", {
  skip_if_not_installed("gt")
  res <- tibble::tibble(
    ref = paste0("tRNA-", 1:5),
    log2FoldChange = rnorm(5),
    pvalue = c(0.01, 0.02, NA, 0.04, NA),
    padj = c(0.05, 0.06, NA, 0.08, NA),
    significant = c(TRUE, FALSE, FALSE, FALSE, FALSE)
  )
  tbl <- tabulate_deseq(res, n = 10)
  expect_equal(nrow(tbl[["_data"]]), 3)
})

test_that("tabulate_deseq works with custom lab_col", {
  skip_if_not_installed("gt")
  res <- tibble::tibble(
    gene = paste0("gene-", 1:5),
    log2FoldChange = rnorm(5),
    pvalue = runif(5, 0, 0.1),
    padj = runif(5, 0, 0.2),
    significant = rep(TRUE, 5)
  )
  tbl <- tabulate_deseq(res, lab_col = "gene")
  expect_s3_class(tbl, "gt_tbl")
  expect_true("gene" %in% names(tbl[["_data"]]))
})
