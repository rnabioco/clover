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

test_that("compute_charging_diffs returns expected columns", {
  config_path <- clover_example("ecoli/config.yaml")
  config <- read_pipeline_config(config_path)
  files <- list_pipeline_files(config, types = "charging")
  paths <- setNames(files$path, files$sample_id)
  charging <- read_charging_multi(paths)
  charging$condition <- ifelse(grepl("ctl", charging$sample_id), "ctl", "inf")

  result <- compute_charging_diffs(
    charging,
    numerator = "inf",
    denominator = "ctl",
    min_count = 50
  )

  expect_s3_class(result, "tbl_df")
  expect_named(
    result,
    c(
      "tRNA",
      "ratio_numerator",
      "ratio_denominator",
      "se_numerator",
      "se_denominator",
      "diff",
      "se_diff"
    )
  )
  expect_true(nrow(result) > 0)
  expect_s3_class(result$tRNA, "factor")
})

test_that("compute_charging_diffs n_top limits rows", {
  config_path <- clover_example("ecoli/config.yaml")
  config <- read_pipeline_config(config_path)
  files <- list_pipeline_files(config, types = "charging")
  paths <- setNames(files$path, files$sample_id)
  charging <- read_charging_multi(paths)
  charging$condition <- ifelse(grepl("ctl", charging$sample_id), "ctl", "inf")

  result <- compute_charging_diffs(
    charging,
    numerator = "inf",
    denominator = "ctl",
    n_top = 5
  )

  expect_lte(nrow(result), 5)
})

test_that("compute_charging_diffs errors when condition column missing", {
  charging <- tibble::tibble(
    tRNA = "tRNA-Ala-AGC-1",
    counts_charged = 100,
    counts_uncharged = 50,
    sample_id = "s1"
  )

  expect_snapshot(
    compute_charging_diffs(charging, numerator = "a", denominator = "b"),
    error = TRUE
  )
})

test_that("pairwise_fisher_exact returns correct structure", {
  mat <- matrix(
    c(1L, 1L, 0L, 0L, 1L, 0L, 1L, 0L, 1L, 0L, 0L, 1L),
    nrow = 6,
    ncol = 2
  )

  result <- pairwise_fisher_exact(mat)

  expect_s3_class(result, "data.frame")
  expect_named(
    result,
    c("pos1", "pos2", "odds_ratio", "log_odds_ratio", "p_value", "total_obs")
  )
  expect_equal(nrow(result), 1)
  expect_equal(result$pos1, 1L)
  expect_equal(result$pos2, 2L)
  expect_equal(result$total_obs, 6L)
})

test_that("pairwise_fisher_exact p-values match fisher.test", {
  mat <- matrix(
    c(
      1L,
      1L,
      1L,
      1L,
      1L,
      0L,
      0L,
      0L,
      0L,
      0L,
      1L,
      1L,
      1L,
      1L,
      0L,
      1L,
      0L,
      0L,
      0L,
      0L
    ),
    nrow = 10,
    ncol = 2
  )

  result <- pairwise_fisher_exact(mat)

  tbl <- table(
    factor(mat[, 1], levels = c(0, 1)),
    factor(mat[, 2], levels = c(0, 1))
  )
  ft <- fisher.test(tbl)

  expect_equal(result$p_value, ft$p.value, tolerance = 1e-6)
  expect_equal(result$odds_ratio, 16)
  expect_equal(result$log_odds_ratio, log(16))
})

test_that("pairwise_fisher_exact handles multiple column pairs", {
  mat <- matrix(
    c(1L, 1L, 0L, 0L, 0L, 0L, 1L, 0L, 1L, 0L, 0L, 0L, 0L, 0L, 0L, 1L, 1L, 0L),
    nrow = 6,
    ncol = 3
  )

  result <- pairwise_fisher_exact(mat)

  expect_equal(nrow(result), 3)
  expect_equal(result$pos1, c(1L, 1L, 2L))
  expect_equal(result$pos2, c(2L, 3L, 3L))
})

test_that("pairwise_fisher_exact skips zero-marginal pairs", {
  mat <- matrix(
    c(0L, 0L, 0L, 0L, 0L, 0L, 1L, 1L, 0L, 0L, 0L, 0L),
    nrow = 6,
    ncol = 2
  )

  result <- pairwise_fisher_exact(mat)

  expect_equal(nrow(result), 0)
})

test_that("pairwise_fisher_exact uses Haldane correction for zero cells", {
  # a=3, b=0, c=2, d=5 -> one cell is zero, Haldane correction applies
  mat <- matrix(
    c(
      1L,
      1L,
      1L,
      0L,
      0L,
      0L,
      0L,
      0L,
      0L,
      0L,
      1L,
      1L,
      1L,
      0L,
      0L,
      1L,
      1L,
      0L,
      0L,
      0L
    ),
    nrow = 10,
    ncol = 2
  )

  result <- pairwise_fisher_exact(mat)

  expected_or <- (3 + 0.5) * (5 + 0.5) / ((0 + 0.5) * (2 + 0.5))
  expect_equal(result$odds_ratio, expected_or)
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
