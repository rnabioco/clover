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

test_that("plot_mod_heatmap works with a single row", {
  df <- data.frame(
    ref = rep("tRNA-Ala", 5),
    sprinzl_label = as.character(1:5),
    value = rnorm(5, sd = 0.1)
  )

  p <- plot_mod_heatmap(df, cluster = TRUE)
  expect_s3_class(p, "ggplot")
})
