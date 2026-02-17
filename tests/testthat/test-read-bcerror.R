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
      "ref", "pos", "cov",
      "a_freq", "t_freq", "g_freq", "c_freq",
      "mis", "ins", "del", "error_rate",
      "mean_qual"
    )
  )

  # ref should be a factor

  expect_true(is.factor(bcerr$ref))
  # pos should be integer
  expect_true(is.integer(bcerr$pos))
  # error_rate, mis, ins, del should be numeric
  expect_true(is.numeric(bcerr$error_rate))
  expect_true(is.numeric(bcerr$mis))

  # 5 charged + 5 uncharged tRNAs in subset data
  expect_equal(length(levels(bcerr$ref)), 10)
})
