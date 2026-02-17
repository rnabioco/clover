test_that("plot_chord_or runs without error", {
  skip_if_not_installed("circlize")

  odds_data <- tibble::tibble(
    pos1 = c(20, 26, 34, 10),
    pos2 = c(34, 44, 58, 25),
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
    pos1 = 20,
    pos2 = 34,
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
    pos1 = c(1, 3),
    pos2 = c(3, 4),
    log_odds_ratio = c(1.2, -0.8),
    p_value = c(0.001, 0.01),
    total_obs = c(100, 200)
  )

  sprinzl <- tibble::tibble(
    trna_id = rep("tRNA-Ala-1", 4),
    seq_index = 1:4,
    sprinzl_label = c("20", "26", "34", "44"),
    global_index = 1:4,
    region = c("D-loop", "anticodon-stem", "anticodon-stem", "variable-region"),
    residue = c("A", "G", "C", "U")
  )

  pdf(nullfile())
  on.exit(dev.off(), add = TRUE)

  result <- plot_chord_or(odds_data, sprinzl_coords = sprinzl)
  expect_null(result)
})

test_that("plot_chord_or converts positions to Sprinzl labels", {
  skip_if_not_installed("circlize")

  odds_data <- tibble::tibble(
    pos1 = c(1, 2),
    pos2 = c(3, 4),
    log_odds_ratio = c(1.2, -0.8),
    p_value = c(0.001, 0.01),
    total_obs = c(100, 200)
  )

  sprinzl <- tibble::tibble(
    trna_id = rep("tRNA-Ala-1", 5),
    seq_index = 1:5,
    sprinzl_label = c("1", "2", "20a", "34", "58"),
    global_index = 1:5,
    region = c(
      "acceptor-stem",
      "acceptor-stem",
      "D-loop",
      "anticodon-stem",
      "T-loop"
    ),
    residue = c("A", "G", "C", "U", "A")
  )

  # Test that .map_to_sprinzl converts correctly
  chord_df <- data.frame(
    from = c("1", "2"),
    to = c("3", "4"),
    value = c(1.2, 0.8),
    stringsAsFactors = FALSE
  )
  mapped <- clover:::.map_to_sprinzl(chord_df, sprinzl)
  expect_equal(mapped$from, c("1", "2"))
  expect_equal(mapped$to, c("20a", "34"))
})

test_that("plot_chord_or shows all Sprinzl positions as sectors", {
  skip_if_not_installed("circlize")

  # Only 2 positions in chord data, but 5 in sprinzl
  odds_data <- tibble::tibble(
    pos1 = 1,
    pos2 = 3,
    log_odds_ratio = 1.5,
    p_value = 0.001,
    total_obs = 100
  )

  sprinzl <- tibble::tibble(
    trna_id = rep("tRNA-Ala-1", 5),
    seq_index = 1:5,
    sprinzl_label = c("1", "2", "3", "4", "5"),
    global_index = 1:5,
    region = rep("acceptor-stem", 5),
    residue = c("A", "G", "C", "U", "A")
  )

  setup <- clover:::.setup_chord_sectors(
    data.frame(from = "1", to = "3", value = 1.5),
    sprinzl
  )
  expect_equal(length(setup$order), 5)
  expect_equal(setup$order, c("1", "2", "3", "4", "5"))
})

test_that(".setup_chord_sectors returns residues", {
  sprinzl <- tibble::tibble(
    trna_id = rep("tRNA-Ala-1", 4),
    seq_index = 1:4,
    sprinzl_label = c("1", "2", "3", "4"),
    global_index = 1:4,
    region = rep("acceptor-stem", 4),
    residue = c("A", "G", "C", "U")
  )

  setup <- clover:::.setup_chord_sectors(
    data.frame(from = "1", to = "3", value = 1.0),
    sprinzl
  )
  expect_equal(setup$residues, c("1" = "A", "2" = "G", "3" = "C", "4" = "U"))
})

test_that("plot_chord_or renders with mods annotation ring", {
  skip_if_not_installed("circlize")

  odds_data <- tibble::tibble(
    pos1 = c(1, 3),
    pos2 = c(3, 4),
    log_odds_ratio = c(1.2, -0.8),
    p_value = c(0.001, 0.01),
    total_obs = c(100, 200)
  )

  sprinzl <- tibble::tibble(
    trna_id = rep("tRNA-Ala-1", 4),
    seq_index = 1:4,
    sprinzl_label = c("20", "26", "34", "44"),
    global_index = 1:4,
    region = c("D-loop", "anticodon-stem", "anticodon-stem", "variable-region"),
    residue = c("A", "G", "C", "U")
  )

  mods <- tibble::tibble(
    ref = rep("tRNA-Ala-1", 2),
    pos = c(1, 3),
    mod_full = c("1-methyladenosine", "pseudouridine"),
    mod1 = c("m1A", "Y")
  )

  pdf(nullfile())
  on.exit(dev.off(), add = TRUE)

  result <- plot_chord_or(
    odds_data,
    sprinzl_coords = sprinzl,
    mods = mods
  )
  expect_null(result)
})

test_that("plot_chord_ror renders with mods annotation ring", {
  skip_if_not_installed("circlize")

  ror_data <- tibble::tibble(
    pos1 = c(1, 3),
    pos2 = c(3, 4),
    or_numerator = c(2.0, -0.5),
    or_denominator = c(1.0, 0.3),
    log_ror = c(1.0, -0.8),
    ror = exp(c(1.0, -0.8))
  )

  sprinzl <- tibble::tibble(
    trna_id = rep("tRNA-Ala-1", 4),
    seq_index = 1:4,
    sprinzl_label = c("20", "26", "34", "44"),
    global_index = 1:4,
    region = c("D-loop", "anticodon-stem", "anticodon-stem", "variable-region"),
    residue = c("A", "G", "C", "U")
  )

  mods <- tibble::tibble(
    ref = rep("tRNA-Ala-1", 2),
    pos = c(1, 3),
    mod_full = c("1-methyladenosine", "pseudouridine"),
    mod1 = c("m1A", "Y")
  )

  pdf(nullfile())
  on.exit(dev.off(), add = TRUE)

  result <- plot_chord_ror(
    ror_data,
    sprinzl_coords = sprinzl,
    mods = mods
  )
  expect_null(result)
})

test_that("compute_ror calculates correct values", {
  odds_data <- tibble::tibble(
    pos1 = rep(20, 4),
    pos2 = rep(34, 4),
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
    pos1 = c(20, 26),
    pos2 = c(34, 44),
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
    pos1 = 20,
    pos2 = 34,
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

test_that(".map_to_sprinzl drops pairs with no mapping", {
  sprinzl <- tibble::tibble(
    trna_id = rep("tRNA-1", 3),
    seq_index = 1:3,
    sprinzl_label = c("1", NA, "3"),
    global_index = 1:3,
    region = rep("D-stem", 3),
    residue = c("A", "G", "C")
  )

  chord_df <- data.frame(
    from = c("1", "1"),
    to = c("2", "3"),
    value = c(1.0, 2.0),
    stringsAsFactors = FALSE
  )

  mapped <- clover:::.map_to_sprinzl(chord_df, sprinzl)
  # Row with to=2 should be dropped (NA sprinzl_label)
  expect_equal(nrow(mapped), 1)
  expect_equal(mapped$from, "1")
  expect_equal(mapped$to, "3")
})
