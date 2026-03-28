test_that("summarize_mod_calls returns expected structure", {
  path <- clover_example("ecoli/mod_calls.tsv.gz")
  res <- summarize_mod_calls(path)

  expect_s3_class(res, "tbl_df")
  expect_named(res, c("ref", "pos", "n_reads", "n_mod", "mod_freq"))
  expect_true(all(res$mod_freq >= 0 & res$mod_freq <= 1))
  expect_true(all(res$n_mod <= res$n_reads))
})

test_that("summarize_mod_calls filters by refs", {
  path <- clover_example("ecoli/mod_calls.tsv.gz")
  res <- summarize_mod_calls(path, refs = "host-tRNA-Glu-TTC-1-1")

  expect_equal(unique(res$ref), "host-tRNA-Glu-TTC-1-1")
})

test_that("summarize_mod_calls filters by min_reads", {
  path <- clover_example("ecoli/mod_calls.tsv.gz")
  res <- summarize_mod_calls(path, min_reads = 20)

  expect_true(all(res$n_reads >= 20))
})

test_that("read_bedmethyl returns expected structure", {
  bed_lines <- c(
    "tRNA-Glu\t10\t11\tm\t50\t+\t10\t11\t255,0,0\t50\t40.0\t20\t30\t0\t0\t0\t0\t0",
    "tRNA-Glu\t20\t21\tm\t5\t+\t20\t21\t255,0,0\t5\t60.0\t3\t2\t0\t0\t0\t0\t0"
  )
  tmp <- withr::local_tempfile(fileext = ".bed")
  writeLines(bed_lines, tmp)

  res <- read_bedmethyl(tmp)

  expect_s3_class(res, "tbl_df")
  expect_named(
    res,
    c("ref", "pos", "strand", "mod_code", "n_valid", "percent_mod",
      "n_mod", "n_canonical")
  )
  expect_equal(nrow(res), 2)
  expect_equal(res$ref, c("tRNA-Glu", "tRNA-Glu"))
  expect_equal(res$pos, c(10, 20))
  expect_equal(res$mod_code, c("m", "m"))
})

test_that("read_bedmethyl filters by min_cov", {
  bed_lines <- c(
    "tRNA-Glu\t10\t11\tm\t50\t+\t10\t11\t255,0,0\t50\t40.0\t20\t30\t0\t0\t0\t0\t0",
    "tRNA-Glu\t20\t21\tm\t5\t+\t20\t21\t255,0,0\t5\t60.0\t3\t2\t0\t0\t0\t0\t0"
  )
  tmp <- withr::local_tempfile(fileext = ".bed")
  writeLines(bed_lines, tmp)

  res <- read_bedmethyl(tmp, min_cov = 10)

  expect_equal(nrow(res), 1)
  expect_equal(res$n_valid, 50L)
})

test_that("summarize_mod_calls output works with compute_bcerror_delta", {
  path <- clover_example("ecoli/mod_calls.tsv.gz")
  mod1 <- summarize_mod_calls(path) |> dplyr::mutate(condition = "wt")
  mod2 <- summarize_mod_calls(path) |> dplyr::mutate(condition = "mut")

  combined <- dplyr::bind_rows(mod1, mod2)
  delta <- compute_bcerror_delta(
    combined,
    delta = wt - mut,
    value_col = "mod_freq"
  )

  expect_s3_class(delta, "tbl_df")
  expect_true("delta" %in% names(delta))
})
