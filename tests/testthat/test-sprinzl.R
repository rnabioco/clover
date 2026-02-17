test_that("read_sprinzl_coords returns expected columns and types", {
  path <- clover_example("sprinzl/sacCer_global_coords.tsv.gz")
  coords <- read_sprinzl_coords(path)

  expect_s3_class(coords, "tbl_df")
  expect_named(
    coords,
    c(
      "trna_id",
      "pos",
      "sprinzl_label",
      "global_index",
      "region",
      "residue"
    )
  )
  expect_type(coords$trna_id, "character")
  expect_type(coords$pos, "double")
  expect_type(coords$sprinzl_label, "character")
  expect_type(coords$global_index, "double")
  expect_type(coords$region, "character")
  expect_type(coords$residue, "character")
})

test_that("read_sprinzl_coords returns reasonable row count for yeast", {
  path <- clover_example("sprinzl/sacCer_global_coords.tsv.gz")
  coords <- read_sprinzl_coords(path)

  # sacCer has ~275 nuclear tRNAs × ~73-93 positions each
  expect_gt(nrow(coords), 10000)
  expect_lt(nrow(coords), 100000)
})

test_that("order_sprinzl_positions orders correctly", {
  labels <- c("20a", "20", "1", "21")
  result <- order_sprinzl_positions(labels)

  expect_s3_class(result, "factor")
  expect_equal(levels(result), c("1", "20", "20a", "21"))
})

test_that("order_sprinzl_positions handles NA and empty labels", {
  labels <- c("1", NA, "", "2", "3")
  result <- order_sprinzl_positions(labels)

  expect_s3_class(result, "factor")
  expect_true(is.na(result[2]))
  # empty string not in levels
  expect_true(is.na(result[3]))
  expect_equal(levels(result), c("1", "2", "3"))
})

test_that("order_sprinzl_positions handles variable loop labels", {
  labels <- c("47:e1", "47:e2", "47", "48")
  result <- order_sprinzl_positions(labels)

  expect_equal(levels(result), c("47", "47:e1", "47:e2", "48"))
})

test_that("trna_regions returns named list of integers", {
  regions <- trna_regions()
  expect_type(regions, "list")
  expect_true(all(vapply(regions, is.integer, logical(1))))
})

test_that("trna_regions contains expected region names", {
  regions <- trna_regions()
  expected <- c(
    "acceptor_stem",
    "d_arm",
    "anticodon_stem",
    "anticodon_loop",
    "variable_loop",
    "t_arm",
    "discriminator",
    "cca"
  )
  expect_named(regions, expected)
})

test_that("trna_regions CCA positions are 74-76", {
  regions <- trna_regions()
  expect_equal(regions$cca, 74L:76L)
})
