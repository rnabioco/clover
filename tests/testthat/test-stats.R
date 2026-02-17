test_that("calc_fold_change computes known values", {
  expect_equal(calc_fold_change(10, 10), 0)
  expect_equal(calc_fold_change(3, 1, pseudocount = 0), log2(3))
})

test_that("calc_fold_change handles zero with pseudocount", {
  result <- calc_fold_change(0, 10, pseudocount = 1)
  expect_equal(result, log2(1 / 11))
})

test_that("calc_fold_change is vectorized", {
  result <- calc_fold_change(c(10, 20), c(5, 10), pseudocount = 0)
  expect_length(result, 2)
  expect_equal(result[1], result[2])
})

test_that("cohens_d returns zero for groups with same mean", {
  set.seed(1)
  x <- rnorm(100, mean = 5, sd = 1)
  y <- rnorm(100, mean = 5, sd = 1)
  expect_equal(cohens_d(x, y), 0, tolerance = 0.3)
})

test_that("cohens_d returns positive value when x > y", {
  set.seed(42)
  x <- rnorm(100, mean = 10, sd = 1)
  y <- rnorm(100, mean = 5, sd = 1)
  expect_gt(cohens_d(x, y), 0)
})

test_that("cohens_d returns known value for unit normal shift", {
  set.seed(1)
  x <- rnorm(10000, mean = 1, sd = 1)
  y <- rnorm(10000, mean = 0, sd = 1)
  expect_equal(cohens_d(x, y), 1, tolerance = 0.1)
})

test_that("propagate_error_ratio computes known case", {
  result <- propagate_error_ratio(10, 5, 1, 0.5)
  expected <- abs(10 / 5) * sqrt((1 / 10)^2 + (0.5 / 5)^2)
  expect_equal(result, expected)
})

test_that("propagate_error_ratio is vectorized", {
  result <- propagate_error_ratio(c(10, 20), c(5, 10), c(1, 2), c(0.5, 1))
  expect_length(result, 2)
})

test_that("propagate_error_diff computes known case", {
  expect_equal(propagate_error_diff(3, 4), 5)
})

test_that("propagate_error_diff with equal errors", {
  expect_equal(propagate_error_diff(1, 1), sqrt(2))
})
