test_that("load_global_coords loads yeast coordinates", {
  coords <- load_global_coords("sacCer")

  expect_s3_class(coords, "tbl_df")
  expect_true(nrow(coords) > 0)
  expect_true(all(c("trna_id", "seq_index", "global_index", "region") %in% names(coords)))

  # Check that global_index is positive integers
  expect_true(all(coords$global_index > 0))

  # Check that regions are present
  expect_true("D-loop" %in% coords$region)
  expect_true("anticodon-loop" %in% coords$region)
})

test_that("add_global_coords joins bcerror with coordinates", {
  bcerr_path <- clover_example("yeast/grande.bcerr.tsv.gz")
  bcerr <- read_bcerror(bcerr_path)

  bcerr_coords <- add_global_coords(bcerr, "sacCer")

  expect_true("global_index" %in% names(bcerr_coords))
  expect_true("sprinzl_label" %in% names(bcerr_coords))
  expect_true("region" %in% names(bcerr_coords))

  # Nuclear tRNAs should have coordinates
  nuc_bcerr <- dplyr::filter(bcerr_coords, grepl("^nuc-", ref))
  expect_true(sum(!is.na(nuc_bcerr$global_index)) > 0)
})

test_that("get_global_labels returns named vector", {
  coords <- load_global_coords("sacCer")
  labels <- get_global_labels(coords)

  expect_type(labels, "character")
  expect_true(length(labels) > 0)
  expect_true(!is.null(names(labels)))
})

test_that("get_region_bounds returns region boundaries", {
  coords <- load_global_coords("sacCer")
  regions <- get_region_bounds(coords)

  expect_s3_class(regions, "tbl_df")
  expect_true(all(c("region", "start", "end") %in% names(regions)))
  expect_true(all(regions$start <= regions$end))
})

test_that("available_organisms returns expected organisms", {
  orgs <- available_organisms()

  expect_true("ecoliK12" %in% orgs)
  expect_true("sacCer" %in% orgs)
  expect_true("hg38" %in% orgs)
})
