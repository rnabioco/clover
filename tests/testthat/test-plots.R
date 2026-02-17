test_that("plot_mod_heatmap returns a ggplot object", {
  df <- tidyr::expand_grid(
    ref = paste0("tRNA-", c("Ala", "Gly", "Ser")),
    sprinzl_label = as.character(1:10)
  )
  df$value <- rnorm(nrow(df), sd = 0.1)

  p <- plot_mod_heatmap(df)
  expect_s3_class(p, "ggplot")
})

test_that("plot_mod_heatmap works with custom column names", {
  df <- tidyr::expand_grid(
    family = paste0("tRNA-", c("Ala", "Gly")),
    sprinzl_label = as.character(1:5)
  )
  df$delta <- rnorm(nrow(df), sd = 0.1)

  p <- plot_mod_heatmap(df, value_col = "delta", ref_col = "family")
  expect_s3_class(p, "ggplot")
})

test_that("plot_mod_heatmap works without clustering", {
  df <- tidyr::expand_grid(
    ref = paste0("tRNA-", c("Ala", "Gly")),
    sprinzl_label = as.character(1:5)
  )
  df$value <- rnorm(nrow(df), sd = 0.1)

  p <- plot_mod_heatmap(df, cluster = FALSE)
  expect_s3_class(p, "ggplot")
})

test_that("plot_volcano returns a ggplot object", {
  skip_if_not_installed("ggrepel")
  res <- tibble::tibble(
    tRNA = paste0("tRNA-", 1:10),
    log2FoldChange = rnorm(10),
    pvalue = c(rep(0.001, 3), rep(0.5, 7)),
    padj = c(rep(0.01, 3), rep(0.8, 7)),
    significant = c(rep(TRUE, 3), rep(FALSE, 7))
  )
  p <- plot_volcano(res)
  expect_s3_class(p, "ggplot")
})

test_that("plot_volcano works with custom lab_col", {
  skip_if_not_installed("ggrepel")
  res <- tibble::tibble(
    gene = paste0("gene-", 1:5),
    log2FoldChange = rnorm(5),
    pvalue = c(0.001, 0.5, 0.5, 0.5, 0.5),
    padj = c(0.01, 0.8, 0.8, 0.8, 0.8),
    significant = c(TRUE, FALSE, FALSE, FALSE, FALSE)
  )
  p <- plot_volcano(res, lab_col = "gene")
  expect_s3_class(p, "ggplot")
})

test_that("plot_volcano works with no significant points", {
  skip_if_not_installed("ggrepel")
  res <- tibble::tibble(
    tRNA = paste0("tRNA-", 1:5),
    log2FoldChange = rnorm(5),
    pvalue = rep(0.5, 5),
    padj = rep(0.8, 5),
    significant = rep(FALSE, 5)
  )
  p <- plot_volcano(res)
  expect_s3_class(p, "ggplot")
})

test_that("plot_mod_heatmap works with a single row", {
  df <- data.frame(
    ref = rep("tRNA-Ala", 5),
    sprinzl_label = as.character(1:5),
    value = rnorm(5, sd = 0.1)
  )

  p <- plot_mod_heatmap(df, cluster = TRUE)
  expect_s3_class(p, "ggplot")
})
