test_that("clean_odds_ratios caps positive infinite values", {
  df <- tibble::tibble(
    ref = rep("tRNA-Ala-AGC-1-1", 3),
    pos1 = c(20, 34, 10),
    pos2 = c(34, 58, 20),
    odds_ratio = c(2.5, Inf, 0.5),
    log_odds_ratio = c(0.92, Inf, -0.69),
    p_value = c(0.001, 0.01, 0.05),
    total_obs = c(200, 150, 100)
  )
  result <- clean_odds_ratios(df)
  expect_true("log_or_clean" %in% names(result))
  expect_true(all(is.finite(result$log_or_clean)))
  expect_gt(
    result$log_or_clean[2],
    max(df$log_odds_ratio[is.finite(df$log_odds_ratio)])
  )
})

test_that("clean_odds_ratios caps negative infinite values", {
  df <- tibble::tibble(
    ref = rep("tRNA-Ala-AGC-1-1", 2),
    pos1 = c(20, 34),
    pos2 = c(34, 58),
    odds_ratio = c(2.5, 0),
    log_odds_ratio = c(0.92, -Inf),
    p_value = c(0.001, 0.01),
    total_obs = c(200, 150)
  )
  result <- clean_odds_ratios(df)
  expect_true(all(is.finite(result$log_or_clean)))
  expect_lt(
    result$log_or_clean[2],
    min(df$log_odds_ratio[is.finite(df$log_odds_ratio)])
  )
})

test_that("clean_odds_ratios preserves finite values", {
  df <- tibble::tibble(
    ref = "tRNA-Ala-AGC-1-1",
    pos1 = 20,
    pos2 = 34,
    odds_ratio = 2.5,
    log_odds_ratio = 0.92,
    p_value = 0.001,
    total_obs = 200
  )
  result <- clean_odds_ratios(df)
  expect_equal(result$log_or_clean, 0.92)
})

test_that("filter_linkages applies all filters", {
  df <- tibble::tibble(
    pos1 = c(20, 34, 10),
    pos2 = c(34, 58, 45),
    odds_ratio = c(4.0, 0.3, 1.1),
    log_odds_ratio = c(1.4, -1.2, 0.1),
    p_adjusted = c(0.001, 0.005, 0.5),
    total_obs = c(200, 150, 50)
  )
  result <- filter_linkages(df)
  expect_equal(nrow(result), 2)
  expect_named(result, c("pos1", "pos2", "value"))
  expect_equal(result$value, c(1.4, -1.2))
})

test_that("filter_linkages respects custom thresholds", {
  df <- tibble::tibble(
    pos1 = c(20, 34),
    pos2 = c(34, 58),
    odds_ratio = c(4.0, 0.3),
    log_odds_ratio = c(1.4, -1.2),
    p_adjusted = c(0.001, 0.005),
    total_obs = c(200, 150)
  )
  result <- filter_linkages(df, max_p = 0.002, min_obs = 100, min_lor = 1.0)
  expect_equal(nrow(result), 1)
  expect_equal(result$pos1, 20)
})

test_that("aggregate_or_isodecoder collapses gene copies", {
  df <- tibble::tibble(
    ref = c("tRNA-Ala-AGC-1-1", "tRNA-Ala-AGC-2-1"),
    pos1 = c(20, 20),
    pos2 = c(34, 34),
    odds_ratio = c(2.0, 4.0),
    log_or_clean = c(0.69, 1.39),
    p_value = c(0.001, 0.01),
    total_obs = c(200, 100)
  )
  result <- aggregate_or_isodecoder(df)
  expect_equal(nrow(result), 1)
  expect_equal(result$isodecoder, "tRNA-Ala-AGC")
  expect_equal(result$n_copies, 2)
  expect_equal(result$mean_or, 3.0)
  expect_equal(result$total_reads, 300)
})

test_that("aggregate_or_isodecoder respects custom pattern", {
  df <- tibble::tibble(
    ref = c("Ala-AGC-copy1", "Ala-AGC-copy2"),
    pos1 = c(20, 20),
    pos2 = c(34, 34),
    odds_ratio = c(2.0, 4.0),
    log_or_clean = c(0.69, 1.39),
    p_value = c(0.001, 0.01),
    total_obs = c(200, 100)
  )
  result <- aggregate_or_isodecoder(df, pattern = "-copy\\d+$")
  expect_equal(result$isodecoder, "Ala-AGC")
})

test_that("compute_ror_isodecoder computes correct ROR", {
  num <- tibble::tibble(
    isodecoder = "tRNA-Ala",
    pos1 = 20,
    pos2 = 34,
    mean_or = 3.0,
    mean_log_or = 1.0,
    sd_log_or = 0.2,
    min_pval = 0.001,
    total_reads = 500,
    n_copies = 2
  )
  den <- tibble::tibble(
    isodecoder = "tRNA-Ala",
    pos1 = 20,
    pos2 = 34,
    mean_or = 1.5,
    mean_log_or = 0.4,
    sd_log_or = 0.15,
    min_pval = 0.01,
    total_reads = 400,
    n_copies = 2
  )
  result <- compute_ror_isodecoder(num, den)
  expect_equal(nrow(result), 1)
  expect_equal(result$ror, 0.6, tolerance = 1e-10)
  expect_equal(result$ror_se, sqrt(0.2^2 + 0.15^2), tolerance = 1e-10)
  expect_true("z_score" %in% names(result))
  expect_true("p_adj" %in% names(result))
  expect_true("significant" %in% names(result))
})

test_that("compute_ror_isodecoder caps extreme values", {
  num <- tibble::tibble(
    isodecoder = "tRNA-Ala",
    pos1 = 20,
    pos2 = 34,
    mean_or = 100,
    mean_log_or = 15,
    sd_log_or = 0.1,
    min_pval = 0.001,
    total_reads = 500,
    n_copies = 2
  )
  den <- tibble::tibble(
    isodecoder = "tRNA-Ala",
    pos1 = 20,
    pos2 = 34,
    mean_or = 0.01,
    mean_log_or = -5,
    sd_log_or = 0.1,
    min_pval = 0.001,
    total_reads = 500,
    n_copies = 2
  )
  result <- compute_ror_isodecoder(num, den, ror_cap = 10)
  expect_equal(result$ror, 10)
})
