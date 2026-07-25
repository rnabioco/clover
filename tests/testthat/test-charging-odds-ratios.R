charging_path <- function() clover_example("ecoli/charging_calls.tsv.gz")
mod_calls_path <- function() clover_example("ecoli/mod_calls.tsv.gz")
mismatch_path <- function() clover_example("ecoli/mismatch_calls.tsv.gz")
bcerror_path <- function() {
  clover_example(
    "ecoli/summary/tables/wt-15-ctl-01/wt-15-ctl-01.bcerror.tsv.gz"
  )
}

test_that("read_charging_calls binarizes at the threshold", {
  result <- read_charging_calls(charging_path())

  expect_s3_class(result, "tbl_df")
  expect_named(result, c("read_id", "ref", "charged"))
  expect_equal(nrow(result), 15)
  expect_equal(sum(result$charged), 6)
  expect_type(result$charged, "integer")
})

test_that("read_charging_calls respects ml_threshold", {
  expect_equal(sum(read_charging_calls(charging_path(), 30)$charged), 15)
  expect_equal(sum(read_charging_calls(charging_path(), 250)$charged), 0)
})

test_that("read_charging_calls errors on a missing column", {
  path <- withr::local_tempfile(fileext = ".tsv")
  readr::write_tsv(tibble::tibble(read_id = "read_1"), path)

  expect_error(read_charging_calls(path), "charging_likelihood")
})

test_that("call_bcerror_sites applies both thresholds", {
  bcerr <- read_bcerror(bcerror_path())
  sites <- call_bcerror_sites(bcerr, min_error = 0.05, min_cov = 10)

  expect_named(sites, c("ref", "pos", "cov", "error_rate"))
  expect_true(all(sites$error_rate >= 0.05))
  expect_true(all(sites$cov >= 10))

  stricter <- call_bcerror_sites(bcerr, min_error = 0.5, min_cov = 10)
  expect_lt(nrow(stricter), nrow(sites))
})

test_that("call_bcerror_sites restricts to refs", {
  bcerr <- read_bcerror(bcerror_path())
  sites <- call_bcerror_sites(
    bcerr,
    min_error = 0.05,
    min_cov = 10,
    refs = "host-tRNA-Glu-TTC-1-1"
  )

  expect_equal(unique(sites$ref), "host-tRNA-Glu-TTC-1-1")
})

test_that("call_bcerror_sites errors on a non-bcerror tibble", {
  expect_error(
    call_bcerror_sites(tibble::tibble(ref = "a", pos = 1L)),
    "missing column"
  )
})

test_that("compute_charging_odds_ratios returns one row per tested site", {
  result <- compute_charging_odds_ratios(
    mod_calls_path(),
    charging_path(),
    min_reads = 5
  )

  expect_s3_class(result, "tbl_df")
  # Positions 34 and 46 carry modification calls; 55 has none, so its
  # modification margin is empty and it is not testable.
  expect_equal(sort(result$pos), c(34L, 46L))
  expect_true(all(result$total_obs == 15))
  expect_true(all(result$p_adjusted >= result$p_value))
})

test_that("compute_charging_odds_ratios computes the 2x2 table correctly", {
  result <- compute_charging_odds_ratios(
    mod_calls_path(),
    charging_path(),
    min_reads = 5
  )
  at34 <- result[result$pos == 34L, ]

  # read_3, read_8 and read_14 are modified at 34 and all three are charged.
  expect_equal(at34$n11, 3)
  expect_equal(at34$n10, 0)
  expect_equal(at34$n01, 3)
  expect_equal(at34$n00, 9)

  # An empty cell triggers the Haldane correction.
  expect_equal(at34$odds_ratio, (3.5 * 9.5) / (0.5 * 3.5), tolerance = 1e-10)
  expect_equal(at34$log_odds_ratio, log(at34$odds_ratio), tolerance = 1e-10)
  expect_equal(at34$mod_freq, 3 / 15, tolerance = 1e-10)
  expect_equal(at34$charged_freq, 6 / 15, tolerance = 1e-10)
})

test_that("compute_charging_odds_ratios p-value matches fisher.test", {
  result <- compute_charging_odds_ratios(
    mod_calls_path(),
    charging_path(),
    min_reads = 5
  )
  at34 <- result[result$pos == 34L, ]

  expected <- stats::fisher.test(matrix(c(3, 3, 0, 9), nrow = 2))$p.value
  expect_equal(at34$p_value, expected, tolerance = 1e-10)
})

test_that("compute_charging_odds_ratios honours min_reads", {
  result <- compute_charging_odds_ratios(
    mod_calls_path(),
    charging_path(),
    min_reads = 100
  )

  expect_equal(nrow(result), 0)
  expect_true(all(
    c("ref", "pos", "odds_ratio", "p_adjusted") %in% names(result)
  ))
})

test_that("compute_charging_odds_ratios accepts tibbles as well as paths", {
  from_paths <- compute_charging_odds_ratios(
    mod_calls_path(),
    charging_path(),
    min_reads = 5
  )
  from_tibbles <- compute_charging_odds_ratios(
    readr::read_tsv(mod_calls_path(), show_col_types = FALSE),
    read_charging_calls(charging_path()),
    min_reads = 5
  )

  expect_equal(from_paths, from_tibbles)
})

test_that("compute_charging_odds_ratios excludes reads lacking a call", {
  calls <- tibble::tibble(
    read_id = c("r1", "r2", "r3", "r4", "r1", "r2"),
    chrom = "tRNA-Ala",
    ref_position = c(10L, 10L, 10L, 10L, 20L, 20L),
    within_alignment = TRUE,
    call_code = c("m", "-", "m", "-", "m", "-")
  )
  charging <- tibble::tibble(
    read_id = c("r1", "r2", "r3", "r4"),
    charged = c(1L, 0L, 1L, 0L)
  )

  result <- compute_charging_odds_ratios(calls, charging, min_reads = 2)

  # Position 20 is covered by two reads only; the other two are missing there
  # rather than unmodified, so they must not inflate total_obs.
  expect_equal(result$total_obs[result$pos == 10L], 4)
  expect_equal(result$total_obs[result$pos == 20L], 2)
})

test_that("compute_charging_odds_ratios collapses repeated calls per read", {
  # modkit emits one row per modification channel, so a read modified at one
  # position can appear several times there. It must still count once.
  calls <- tibble::tibble(
    read_id = c("r1", "r1", "r2", "r3", "r4"),
    chrom = "tRNA-Ala",
    ref_position = 10L,
    within_alignment = TRUE,
    call_code = c("m", "a", "-", "m", "-")
  )
  charging <- tibble::tibble(
    read_id = c("r1", "r2", "r3", "r4"),
    charged = c(1L, 0L, 1L, 0L)
  )

  result <- compute_charging_odds_ratios(calls, charging, min_reads = 2)

  expect_equal(result$total_obs, 4)
  expect_equal(result$n11, 2)
  expect_equal(result$n00, 2)
})

test_that("compute_charging_odds_ratios dedupe = FALSE skips the collapse", {
  # One row per read and position, as the pipeline's mismatch calls emit.
  calls <- tibble::tibble(
    read_id = c("r1", "r2", "r3", "r4"),
    chrom = "tRNA-Ala",
    ref_position = 10L,
    within_alignment = TRUE,
    call_code = c("X", "-", "X", "-")
  )
  charging <- tibble::tibble(
    read_id = c("r1", "r2", "r3", "r4"),
    charged = c(1L, 0L, 1L, 0L)
  )

  expect_equal(
    compute_charging_odds_ratios(calls, charging, min_reads = 2, dedupe = FALSE),
    compute_charging_odds_ratios(calls, charging, min_reads = 2, dedupe = TRUE)
  )
})

test_that("compute_charging_odds_ratios prunes thin margins", {
  # 20 reads, but only one of them is modified: the odds ratio is defined and
  # the site is testable in principle, yet it carries no useful power.
  calls <- tibble::tibble(
    read_id = paste0("r", 1:20),
    chrom = "tRNA-Ala",
    ref_position = 10L,
    within_alignment = TRUE,
    call_code = c("m", rep("-", 19))
  )
  charging <- tibble::tibble(
    read_id = paste0("r", 1:20),
    charged = rep(c(1L, 0L), each = 10)
  )

  expect_equal(nrow(compute_charging_odds_ratios(calls, charging)), 1)
  expect_equal(
    nrow(compute_charging_odds_ratios(calls, charging, min_margin = 5)),
    0
  )
})

test_that("compute_charging_odds_ratios prunes sites that cannot reach max_p", {
  # Six reads split 3/3 on charging with one modified read: the most extreme
  # outcome available still cannot clear 0.05.
  calls <- tibble::tibble(
    read_id = paste0("r", 1:6),
    chrom = "tRNA-Ala",
    ref_position = 10L,
    within_alignment = TRUE,
    call_code = c("m", rep("-", 5))
  )
  charging <- tibble::tibble(
    read_id = paste0("r", 1:6),
    charged = rep(c(1L, 0L), each = 3)
  )

  unpruned <- compute_charging_odds_ratios(calls, charging, min_reads = 5)
  expect_equal(nrow(unpruned), 1)
  expect_gt(unpruned$p_value, 0.05)

  expect_equal(
    nrow(compute_charging_odds_ratios(
      calls, charging, min_reads = 5, max_p = 0.05
    )),
    0
  )
})

test_that("compute_charging_odds_ratios max_p keeps sites that can reach it", {
  # A site with the same read count but a balanced modification margin can
  # reach significance, so it must survive the same filter.
  calls <- tibble::tibble(
    read_id = paste0("r", 1:20),
    chrom = "tRNA-Ala",
    ref_position = 10L,
    within_alignment = TRUE,
    call_code = rep(c("m", "-"), each = 10)
  )
  charging <- tibble::tibble(
    read_id = paste0("r", 1:20),
    charged = rep(c(1L, 0L), each = 10)
  )

  result <- compute_charging_odds_ratios(
    calls, charging, min_reads = 5, max_p = 0.05
  )

  expect_equal(nrow(result), 1)
  expect_lt(result$p_value, 0.05)
})

test_that("compute_charging_odds_ratios pruning does not alter kept results", {
  calls <- readr::read_tsv(mod_calls_path(), show_col_types = FALSE)
  charging <- read_charging_calls(charging_path())

  full <- compute_charging_odds_ratios(calls, charging, min_reads = 5)
  pruned <- compute_charging_odds_ratios(
    calls, charging, min_reads = 5, min_margin = 2
  )

  common <- dplyr::semi_join(full, pruned, by = c("ref", "pos"))
  expect_equal(common$p_value, pruned$p_value)
  expect_equal(common$odds_ratio, pruned$odds_ratio)
})

test_that("compute_charging_odds_ratios drops reads with no charging call", {
  calls <- tibble::tibble(
    read_id = c("r1", "r2", "r3"),
    chrom = "tRNA-Ala",
    ref_position = 10L,
    within_alignment = TRUE,
    call_code = c("m", "-", "m")
  )
  charging <- tibble::tibble(read_id = c("r1", "r2"), charged = c(1L, 0L))

  result <- compute_charging_odds_ratios(calls, charging, min_reads = 2)

  expect_equal(result$total_obs, 2)
})

test_that("compute_charging_odds_ratios restricts to bcerror sites", {
  sites <- call_bcerror_sites(
    read_bcerror(bcerror_path()),
    min_error = 0.05,
    min_cov = 10,
    refs = "host-tRNA-Glu-TTC-1-1"
  )
  result <- compute_charging_odds_ratios(
    mismatch_path(),
    charging_path(),
    sites = sites,
    min_reads = 5
  )

  expect_true(all(result$pos %in% sites$pos))
  expect_setequal(result$pos, c(34L, 55L, 58L))

  # Position 58 mismatches track charging closely; 55 runs the other way.
  expect_gt(result$log_odds_ratio[result$pos == 58L], 0)
  expect_lt(result$log_odds_ratio[result$pos == 55L], 0)
})

test_that("compute_charging_odds_ratios drops sites outside the site list", {
  sites <- tibble::tibble(ref = "host-tRNA-Glu-TTC-1-1", pos = 34L)
  result <- compute_charging_odds_ratios(
    mod_calls_path(),
    charging_path(),
    sites = sites,
    min_reads = 5
  )

  expect_equal(result$pos, 34L)
})

test_that("compute_charging_odds_ratios warns when no site overlaps", {
  sites <- tibble::tibble(ref = "host-tRNA-Glu-TTC-1-1", pos = 9999L)

  expect_warning(
    result <- compute_charging_odds_ratios(
      mod_calls_path(),
      charging_path(),
      sites = sites
    ),
    "No calls overlap"
  )
  expect_equal(nrow(result), 0)
})

test_that("compute_charging_odds_ratios restricts to refs", {
  result <- compute_charging_odds_ratios(
    mod_calls_path(),
    charging_path(),
    refs = "no-such-tRNA"
  )

  expect_equal(nrow(result), 0)
})

test_that("compute_charging_odds_ratios errors on malformed input", {
  charging <- tibble::tibble(read_id = "r1", charged = 1L)

  expect_error(
    compute_charging_odds_ratios(
      tibble::tibble(read_id = "r1", chrom = "a", ref_position = 1L),
      charging
    ),
    "call_code"
  )
  expect_error(
    compute_charging_odds_ratios(mod_calls_path(), tibble::tibble(x = 1)),
    "charging_likelihood"
  )
})

test_that("compute_odds_ratios restricts pairs to supplied sites", {
  sites <- tibble::tibble(
    ref = "host-tRNA-Glu-TTC-1-1",
    pos = c(34L, 46L)
  )
  result <- compute_odds_ratios(mod_calls_path(), min_reads = 5, sites = sites)

  expect_equal(nrow(result), 1)
  expect_setequal(c(result$pos1, result$pos2), c("34", "46"))
})
