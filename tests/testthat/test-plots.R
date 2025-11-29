# Tests for plotting functions

# Test classify_trna_type helper ----
test_that("classify_trna_type identifies Type II correctly", {
  type2_ids <- c("tRNA-Leu-CAA-1-1", "tRNA-Ser-GCT-1-1", "tRNA-Tyr-GTA-1-1")
  type1_ids <- c("tRNA-Ala-GGC-1-1", "tRNA-Arg-ACG-1-1", "tRNA-Gly-GCC-1-1")

  expect_equal(classify_trna_type(type2_ids), rep("Type II", 3))
  expect_equal(classify_trna_type(type1_ids), rep("Type I", 3))
})

test_that("classify_trna_type treats SeC as Type I (excluded from coordinates)", {
  # SeC tRNAs are excluded from coordinate system due to structural incompatibility
  # classify_trna_type returns "Type I" (default) but SeC won't match any group
  sec_id <- "tRNA-SeC-TCA-1-1"
  expect_equal(classify_trna_type(sec_id), "Type I")
})

# Test plot_bcerror_heatmap split behavior ----
test_that("plot_bcerror_heatmap returns patchwork object when split", {
  bcerr_coords <- read_bcerror(clover_example("yeast/grande.bcerr.tsv.gz")) |>
    add_global_coords("sacCer")

  p <- plot_bcerror_heatmap(bcerr_coords, split_by_group = TRUE)

  # Should return patchwork object
  expect_s3_class(p, "patchwork")
})

test_that("plot_bcerror_heatmap returns single ggplot when not split", {
  bcerr_coords <- read_bcerror(clover_example("yeast/grande.bcerr.tsv.gz")) |>
    add_global_coords("sacCer")

  p <- plot_bcerror_heatmap(bcerr_coords, split_by_group = FALSE)

  # Should return single ggplot
  expect_s3_class(p, "gg")
  expect_s3_class(p, "ggplot")
})

test_that("plot_bcerror_heatmap shows Type I and Type II separately", {
  bcerr_coords <- read_bcerror(clover_example("yeast/grande.bcerr.tsv.gz")) |>
    add_global_coords("sacCer") |>
    dplyr::filter(grepl("^nuc-", ref))

  # Add type classification to check
  bcerr_with_type <- bcerr_coords |>
    dplyr::mutate(trna_type = classify_trna_type(ref))

  # Verify we have both types in yeast data
  expect_true("Type I" %in% bcerr_with_type$trna_type)
  expect_true("Type II" %in% bcerr_with_type$trna_type)

  # Plot should work
  p <- plot_bcerror_heatmap(bcerr_coords, split_by_group = TRUE)
  expect_s3_class(p, "patchwork")
})

test_that("plot_bcerror_heatmap requires global coordinates", {
  bcerr <- read_bcerror(clover_example("yeast/grande.bcerr.tsv.gz"))
  # Don't add global coordinates

  expect_error(
    plot_bcerror_heatmap(bcerr),
    "global coordinates"
  )
})

# Test axis labeling ----
test_that("plot_bcerror_heatmap shows all Sprinzl labels", {
  bcerr_coords <- read_bcerror(clover_example("yeast/grande.bcerr.tsv.gz")) |>
    add_global_coords("sacCer") |>
    dplyr::filter(grepl("^nuc-", ref))

  # Get single plot to test labels
  p <- plot_bcerror_heatmap(bcerr_coords, split_by_group = FALSE)

  # Build plot to access x-axis
  built <- ggplot2::ggplot_build(p)
  x_breaks <- built$layout$panel_params[[1]]$x$breaks

  # Should have many labels (Type I tRNAs have ~76 positions)
  # With all labels shown, expect > 50 positions
  expect_gt(length(x_breaks), 50)
})

test_that("Internal function converts -1 to NA in labels", {
  bcerr_coords <- read_bcerror(clover_example("yeast/grande.bcerr.tsv.gz")) |>
    add_global_coords("sacCer") |>
    dplyr::filter(!is.na(global_index)) |>
    dplyr::filter(grepl("^nuc-", ref))

  p <- plot_bcerror_heatmap(bcerr_coords, split_by_group = FALSE)

  # Plot should build without errors (implementation converts "-1" to "NA")
  built <- ggplot2::ggplot_build(p)
  expect_s3_class(built, "ggplot_built")
})
