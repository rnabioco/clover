test_that("plot_chord_or runs without error", {
  skip_if_not_installed("circlize")

  odds_data <- tibble::tibble(
    pos1 = c("20", "26", "34", "10"),
    pos2 = c("34", "44", "58", "25"),
    log_odds_ratio = c(1.2, -0.8, 0.6, -1.5),
    p_value = c(0.001, 0.01, 0.04, 0.001),
    total_obs = c(100, 200, 150, 80)
  )

  pdf(nullfile())
  on.exit(dev.off(), add = TRUE)

  result <- plot_chord_or(odds_data, or_cutoff = 0.5)
  expect_null(result)
})

test_that("plot_chord_or returns NULL when no significant pairs", {
  skip_if_not_installed("circlize")

  odds_data <- tibble::tibble(
    pos1 = "20",
    pos2 = "34",
    log_odds_ratio = 0.1,
    p_value = 0.5,
    total_obs = 100
  )

  pdf(nullfile())
  on.exit(dev.off(), add = TRUE)

  expect_message(
    result <- plot_chord_or(odds_data),
    "No significant pairs"
  )
  expect_null(result)
})

test_that("plot_chord_or works with sprinzl_coords", {
  skip_if_not_installed("circlize")

  odds_data <- tibble::tibble(
    pos1 = c("20", "26"),
    pos2 = c("34", "44"),
    log_odds_ratio = c(1.2, -0.8),
    p_value = c(0.001, 0.01),
    total_obs = c(100, 200)
  )

  sprinzl <- tibble::tibble(
    trna_id = rep("tRNA-Ala-1", 4),
    seq_index = 1:4,
    sprinzl_label = c("20", "26", "34", "44"),
    global_index = 1:4,
    region = c("D-loop", "AC-stem", "AC-stem", "variable-loop"),
    residue = c("A", "G", "C", "U")
  )

  pdf(nullfile())
  on.exit(dev.off(), add = TRUE)

  result <- plot_chord_or(odds_data, sprinzl_coords = sprinzl)
  expect_null(result)
})

test_that("compute_ror calculates correct values", {
  odds_data <- tibble::tibble(
    pos1 = rep("20", 4),
    pos2 = rep("34", 4),
    log_odds_ratio = c(1.0, 1.2, 2.0, 2.4),
    total_obs = rep(100, 4),
    sample_id = c("wt1", "wt2", "mut1", "mut2"),
    condition = rep(c("wt", "mut"), each = 2)
  )

  result <- compute_ror(
    odds_data,
    numerator = "mut",
    denominator = "wt"
  )

  expect_s3_class(result, "tbl_df")
  expect_named(
    result,
    c("pos1", "pos2", "or_numerator", "or_denominator", "log_ror", "ror")
  )
  expect_equal(nrow(result), 1)
  # mean(2.0, 2.4) - mean(1.0, 1.2) = 2.2 - 1.1 = 1.1
  expect_equal(result$or_numerator, 2.2)
  expect_equal(result$or_denominator, 1.1)
  expect_equal(result$log_ror, 1.1, tolerance = 1e-10)
  expect_equal(result$ror, exp(1.1), tolerance = 1e-10)
})

test_that("plot_chord_ror runs without error", {
  skip_if_not_installed("circlize")

  ror_data <- tibble::tibble(
    pos1 = c("20", "26"),
    pos2 = c("34", "44"),
    or_numerator = c(2.0, -0.5),
    or_denominator = c(1.0, 0.3),
    log_ror = c(1.0, -0.8),
    ror = exp(c(1.0, -0.8))
  )

  pdf(nullfile())
  on.exit(dev.off(), add = TRUE)

  result <- plot_chord_ror(ror_data, ror_cutoff = 0.5)
  expect_null(result)
})

test_that("plot_chord_ror returns NULL when no pairs exceed cutoff", {
  skip_if_not_installed("circlize")

  ror_data <- tibble::tibble(
    pos1 = "20",
    pos2 = "34",
    log_ror = 0.1,
    ror = exp(0.1)
  )

  pdf(nullfile())
  on.exit(dev.off(), add = TRUE)

  expect_message(
    result <- plot_chord_ror(ror_data, ror_cutoff = 0.5),
    "No pairs exceed"
  )
  expect_null(result)
})
