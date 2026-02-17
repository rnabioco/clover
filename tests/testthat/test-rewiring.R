test_that("prepare_rewiring_matrix returns expected matrix", {
  df <- tibble::tibble(
    isodecoder = rep(c("tRNA-Ala", "tRNA-Gly"), each = 2),
    pos1 = c(20, 34, 20, 34),
    pos2 = c(34, 58, 34, 58),
    ror = c(1.5, -0.8, 0.3, 2.1),
    significant = c(TRUE, TRUE, TRUE, TRUE)
  )
  mat <- prepare_rewiring_matrix(df)
  expect_true(is.matrix(mat))
  expect_equal(nrow(mat), 2)
  expect_equal(ncol(mat), 2)
  expect_equal(sort(rownames(mat)), c("tRNA-Ala", "tRNA-Gly"))
})

test_that("prepare_rewiring_matrix filters by significance", {
  df <- tibble::tibble(
    isodecoder = rep("tRNA-Ala", 2),
    pos1 = c(20, 34),
    pos2 = c(34, 58),
    ror = c(1.5, -0.8),
    significant = c(TRUE, FALSE)
  )
  mat <- prepare_rewiring_matrix(df, sig_only = TRUE)
  expect_equal(ncol(mat), 1)

  mat_all <- prepare_rewiring_matrix(df, sig_only = FALSE)
  expect_equal(ncol(mat_all), 2)
})

test_that("prepare_rewiring_matrix caps extreme values", {
  df <- tibble::tibble(
    isodecoder = "tRNA-Ala",
    pos1 = 20,
    pos2 = 34,
    ror = 15,
    significant = TRUE
  )
  mat <- prepare_rewiring_matrix(df, value_cap = 10)
  expect_equal(mat[1, 1], 10)
})

test_that("calculate_rewiring_scores returns correct columns", {
  mat <- matrix(
    c(1.5, -0.8, 0.3, 2.1),
    nrow = 2,
    dimnames = list(c("tRNA-Ala", "tRNA-Gly"), c("20_vs_34", "34_vs_58"))
  )
  scores <- calculate_rewiring_scores(mat)
  expect_s3_class(scores, "tbl_df")
  expect_named(
    scores,
    c(
      "isodecoder",
      "euclidean_magnitude",
      "mean_abs_change",
      "max_abs_change",
      "n_nonzero"
    )
  )
  expect_equal(nrow(scores), 2)
})

test_that("calculate_rewiring_scores computes euclidean magnitude", {
  mat <- matrix(
    c(3, 4),
    nrow = 1,
    dimnames = list("tRNA-Ala", c("a", "b"))
  )
  scores <- calculate_rewiring_scores(mat)
  expect_equal(scores$euclidean_magnitude, 5, ignore_attr = TRUE)
})

test_that("calculate_rewiring_scores sorts by magnitude", {
  mat <- matrix(
    c(1, 10, 1, 1),
    nrow = 2,
    dimnames = list(c("small", "big"), c("a", "b"))
  )
  scores <- calculate_rewiring_scores(mat)
  expect_equal(scores$isodecoder[1], "big")
})

test_that("perform_pcoa returns expected structure", {
  mat <- matrix(
    c(1.5, -0.8, 0.3, 2.1, 0.5, -1.2),
    nrow = 3,
    dimnames = list(
      c("tRNA-Ala", "tRNA-Gly", "tRNA-Ser"),
      c("20_vs_34", "34_vs_58")
    )
  )
  result <- perform_pcoa(mat)
  expect_type(result, "list")
  expect_named(result, c("coordinates", "variance_explained", "eigenvalues"))
  expect_s3_class(result$coordinates, "tbl_df")
  expect_true("isodecoder" %in% names(result$coordinates))
  expect_true("PC1" %in% names(result$coordinates))
  expect_true("PC2" %in% names(result$coordinates))
  expect_length(result$variance_explained, 2)
})

test_that("perform_pcoa respects k parameter", {
  mat <- matrix(
    rnorm(12),
    nrow = 4,
    dimnames = list(paste0("t", 1:4), paste0("c", 1:3))
  )
  result <- perform_pcoa(mat, k = 3)
  expect_true("PC3" %in% names(result$coordinates))
  expect_length(result$variance_explained, 3)
})
