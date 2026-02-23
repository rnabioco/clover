test_that("read_bcerror works with new pipeline format", {
  bcerr_path <- clover_example(
    "ecoli/summary/tables/wt-15-ctl-01/wt-15-ctl-01.bcerror.tsv.gz"
  )
  bcerr <- read_bcerror(bcerr_path)

  expect_s3_class(bcerr, "tbl_df")
  expect_true(nrow(bcerr) > 0)
  expect_named(
    bcerr,
    c(
      "ref",
      "pos",
      "cov",
      "a_freq",
      "t_freq",
      "g_freq",
      "c_freq",
      "mis",
      "ins",
      "del",
      "error_rate",
      "mean_qual"
    )
  )

  # ref should be character
  expect_type(bcerr$ref, "character")
  # pos should be integer
  expect_true(is.integer(bcerr$pos))
  # error_rate, mis, ins, del should be numeric
  expect_true(is.numeric(bcerr$error_rate))
  expect_true(is.numeric(bcerr$mis))

  # 5 charged + 5 uncharged tRNAs in subset data
  expect_equal(length(unique(bcerr$ref)), 10)
})

test_that("compute_bcerror_delta computes lhs - rhs", {
  df <- tidyr::expand_grid(
    ref = c("tRNA-Ala", "tRNA-Gly"),
    pos = 1:3,
    condition = c("wt", "mut")
  )
  df$mean_error <- c(
    0.1, 0.2, 0.3, 0.05, 0.15, 0.25,
    0.4, 0.5, 0.6, 0.35, 0.45, 0.55
  )

  result <- compute_bcerror_delta(df, delta = wt - mut)

  expect_s3_class(result, "tbl_df")
  expect_named(result, c("ref", "pos", "wt", "mut", "delta"))
  expect_equal(nrow(result), 6)
  expect_equal(result$delta, result$wt - result$mut)
})

test_that("compute_bcerror_delta errors on invalid delta expression", {
  df <- tibble::tibble(
    ref = "tRNA-Ala", pos = 1L, condition = "wt", mean_error = 0.1
  )
  expect_snapshot(compute_bcerror_delta(df, delta = wt + mut), error = TRUE)
})

test_that("compute_bcerror_delta errors on missing condition level", {
  df <- tibble::tibble(
    ref = "tRNA-Ala", pos = 1L, condition = "wt", mean_error = 0.1
  )
  expect_snapshot(compute_bcerror_delta(df, delta = wt - missing), error = TRUE)
})

test_that("compute_bcerror_delta respects custom column names", {
  df <- tibble::tibble(
    ref = rep("tRNA-Ala", 2),
    pos = c(1L, 1L),
    group = c("a", "b"),
    err = c(0.3, 0.1)
  )
  result <- compute_bcerror_delta(
    df,
    delta = a - b,
    value_col = "err",
    condition_col = "group"
  )
  expect_equal(result$delta, 0.2)
})
