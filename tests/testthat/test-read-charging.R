test_that("read_charging reads a TSV file", {
  tmp <- withr::local_tempfile(fileext = ".tsv")
  writeLines(
    c(
      "tRNA\tcounts_charged\tcounts_uncharged\tcpm_charged\tcpm_uncharged",
      "tRNA-Ala-AGC-1\t100\t50\t666.7\t333.3",
      "tRNA-Gly-GCC-1\t200\t80\t714.3\t285.7"
    ),
    tmp
  )

  result <- read_charging(tmp)

  expect_s3_class(result, "tbl_df")
  expect_equal(nrow(result), 2)
  expect_true(all(
    c("tRNA", "counts_charged", "counts_uncharged") %in% names(result)
  ))
  expect_equal(result$counts_charged, c(100, 200))
})

test_that("read_odds_ratios reads a TSV file", {
  tmp <- withr::local_tempfile(fileext = ".tsv")
  writeLines(
    c(
      "pos1\tpos2\todds_ratio\tlog_odds_ratio\tp_value\ttotal_obs",
      "20\t34\t2.5\t0.916\t0.01\t100",
      "26\t44\t0.5\t-0.693\t0.03\t200"
    ),
    tmp
  )

  result <- read_odds_ratios(tmp)

  expect_s3_class(result, "tbl_df")
  expect_equal(nrow(result), 2)
  expect_true(all(
    c("pos1", "pos2", "odds_ratio", "log_odds_ratio") %in% names(result)
  ))
})

test_that("read_charging_multi combines samples", {
  tmp1 <- withr::local_tempfile(fileext = ".tsv")
  tmp2 <- withr::local_tempfile(fileext = ".tsv")

  header <- "tRNA\tcounts_charged\tcounts_uncharged"
  writeLines(c(header, "tRNA-Ala\t100\t50"), tmp1)
  writeLines(c(header, "tRNA-Ala\t120\t60"), tmp2)

  paths <- c(sample1 = tmp1, sample2 = tmp2)
  result <- read_charging_multi(paths)

  expect_s3_class(result, "tbl_df")
  expect_equal(nrow(result), 2)
  expect_true("sample_id" %in% names(result))
  expect_equal(sort(unique(result$sample_id)), c("sample1", "sample2"))
})

test_that("read_charging_multi errors on unnamed paths", {
  expect_error(
    read_charging_multi(c("/fake/path1", "/fake/path2")),
    "named character vector"
  )
})

test_that("read_odds_ratios_multi combines samples", {
  tmp1 <- withr::local_tempfile(fileext = ".tsv")
  tmp2 <- withr::local_tempfile(fileext = ".tsv")

  header <- "pos1\tpos2\todds_ratio\tlog_odds_ratio\tp_value\ttotal_obs"
  writeLines(c(header, "20\t34\t2.5\t0.916\t0.01\t100"), tmp1)
  writeLines(c(header, "20\t34\t3.0\t1.099\t0.005\t150"), tmp2)

  paths <- c(wt = tmp1, mut = tmp2)
  result <- read_odds_ratios_multi(paths)

  expect_s3_class(result, "tbl_df")
  expect_equal(nrow(result), 2)
  expect_true("sample_id" %in% names(result))
  expect_equal(sort(unique(result$sample_id)), c("mut", "wt"))
})

test_that("read_odds_ratios_multi errors on unnamed paths", {
  expect_error(
    read_odds_ratios_multi(c("/fake/path1")),
    "named character vector"
  )
})

test_that("read_charging works with ecoli test data", {
  path <- clover_example(
    "ecoli/summary/tables/wt-15-ctl-01/wt-15-ctl-01.charging.cpm.tsv.gz"
  )
  result <- read_charging(path)

  expect_s3_class(result, "tbl_df")
  expect_true(nrow(result) > 0)
  expect_true(all(
    c("tRNA", "counts_charged", "counts_uncharged") %in% names(result)
  ))
})

test_that("read_odds_ratios works with ecoli test data", {
  path <- clover_example(
    "ecoli/summary/tables/wt-15-ctl-01/wt-15-ctl-01.odds_ratios.tsv.gz"
  )
  result <- read_odds_ratios(path)

  expect_s3_class(result, "tbl_df")
  expect_true(nrow(result) > 0)
  expect_true(all(
    c("ref", "pos1", "pos2", "odds_ratio", "log_odds_ratio") %in%
      names(result)
  ))
})
