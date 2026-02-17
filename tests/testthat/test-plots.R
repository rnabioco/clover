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

test_that("plot_abundance_charging returns a ggplot object", {
  skip_if_not_installed("ggrepel")
  deseq_res <- tibble::tibble(
    tRNA = paste0("tRNA-", 1:6),
    log2FoldChange = c(1, -1, 0.5, -0.5, 2, -2),
    padj = c(0.01, 0.02, 0.5, 0.6, 0.001, 0.003)
  )
  charging_diffs <- tibble::tibble(
    tRNA = paste0("tRNA-", 1:6),
    diff = c(0.1, -0.1, 0.05, -0.05, -0.2, 0.15),
    se_diff = rep(0.03, 6)
  )
  p <- plot_abundance_charging(deseq_res, charging_diffs)
  expect_s3_class(p, "ggplot")
})

test_that("plot_abundance_charging handles no significant points", {
  skip_if_not_installed("ggrepel")
  deseq_res <- tibble::tibble(
    tRNA = paste0("tRNA-", 1:5),
    log2FoldChange = rnorm(5),
    padj = rep(0.8, 5)
  )
  charging_diffs <- tibble::tibble(
    tRNA = paste0("tRNA-", 1:5),
    diff = rnorm(5, sd = 0.1),
    se_diff = rep(0.03, 5)
  )
  p <- plot_abundance_charging(deseq_res, charging_diffs)
  expect_s3_class(p, "ggplot")
})

test_that("plot_abundance_charging respects custom padj_cutoff", {
  skip_if_not_installed("ggrepel")
  deseq_res <- tibble::tibble(
    tRNA = paste0("tRNA-", 1:5),
    log2FoldChange = c(1, -1, 0.5, -0.5, 2),
    padj = c(0.005, 0.02, 0.05, 0.1, 0.001)
  )
  charging_diffs <- tibble::tibble(
    tRNA = paste0("tRNA-", 1:5),
    diff = c(0.1, -0.1, 0.05, -0.05, -0.2),
    se_diff = rep(0.03, 5)
  )
  p_strict <- plot_abundance_charging(
    deseq_res,
    charging_diffs,
    padj_cutoff = 0.01
  )
  sig_strict <- sum(p_strict$data$significant)
  p_loose <- plot_abundance_charging(
    deseq_res,
    charging_diffs,
    padj_cutoff = 0.1
  )
  sig_loose <- sum(p_loose$data$significant)
  expect_gt(sig_loose, sig_strict)
})

test_that("plot_charging_diffs returns a ggplot object", {
  df <- tibble::tibble(
    tRNA = forcats::fct_inorder(paste0("tRNA-", 1:5)),
    diff = c(-0.1, -0.05, 0.02, 0.08, 0.15),
    se_diff = rep(0.03, 5)
  )
  p <- plot_charging_diffs(df)
  expect_s3_class(p, "ggplot")
})

test_that("plot_charging_diffs respects point_size", {
  df <- tibble::tibble(
    tRNA = forcats::fct_inorder(paste0("tRNA-", 1:3)),
    diff = c(-0.1, 0.0, 0.1),
    se_diff = rep(0.02, 3)
  )
  p <- plot_charging_diffs(df, point_size = 4)
  expect_s3_class(p, "ggplot")
})

test_that("plot_bcerror_profile returns a ggplot object", {
  df <- tidyr::expand_grid(
    ref = c("tRNA-Ala", "tRNA-Gly"),
    pos = 1:20,
    condition = c("ctl", "inf")
  )
  df$mean_error <- runif(nrow(df), 0, 0.3)
  p <- plot_bcerror_profile(df)
  expect_s3_class(p, "ggplot")
})

test_that("plot_bcerror_profile filters by refs", {
  df <- tidyr::expand_grid(
    ref = c("tRNA-Ala", "tRNA-Gly", "tRNA-Ser"),
    pos = 1:10,
    condition = c("ctl", "inf")
  )
  df$mean_error <- runif(nrow(df), 0, 0.3)
  p <- plot_bcerror_profile(df, refs = c("tRNA-Ala", "tRNA-Gly"))
  expect_s3_class(p, "ggplot")
  expect_equal(length(unique(p$data$ref)), 2)
})

test_that("plot_bcerror_profile overlays modification positions", {
  df <- tidyr::expand_grid(
    ref = "tRNA-Ala",
    pos = 1:20,
    condition = c("ctl", "inf")
  )
  df$mean_error <- runif(nrow(df), 0, 0.3)
  mods <- tibble::tibble(ref = "tRNA-Ala", pos = c(5, 10, 15))
  p <- plot_bcerror_profile(df, mods = mods)
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

test_that("plot_pcoa_rewiring returns a ggplot object", {
  skip_if_not_installed("ggrepel")

  mat <- matrix(
    c(1.5, -0.8, 0.3, 2.1, 0.5, -1.2),
    nrow = 3,
    dimnames = list(
      c("tRNA-Ala", "tRNA-Gly", "tRNA-Ser"),
      c("20_vs_34", "34_vs_58")
    )
  )
  scores <- calculate_rewiring_scores(mat)
  pcoa <- perform_pcoa(mat)
  p <- plot_pcoa_rewiring(pcoa, scores)
  expect_s3_class(p, "ggplot")
})
