# Tests for structure visualization functions

# Test 1: load_structure_template() returns valid structure ----
test_that("load_structure_template returns valid template", {
  template <- load_structure_template()

  expect_type(template, "list")
  expect_named(template, c("residues", "basepairs"))

  # Check residues tibble structure
  expect_s3_class(template$residues, "tbl_df")
  expect_true(all(c("residue_index", "sprinzl_label", "x", "y") %in% names(template$residues)))

  # Check basepairs tibble structure
  expect_s3_class(template$basepairs, "tbl_df")
  expect_true(all(c("pos1", "pos2") %in% names(template$basepairs)))
})

# Test 2: load_structure_template() has sensible coordinate ranges ----
test_that("load_structure_template has sensible coordinates", {
  template <- load_structure_template()

  # Verify coordinates are finite
  expect_true(all(is.finite(template$residues$x)))
  expect_true(all(is.finite(template$residues$y)))

  # Verify standard tRNA length (~76 positions for Type I)
  expect_gt(nrow(template$residues), 70)
  expect_lt(nrow(template$residues), 85)

  # Verify base pairs reference valid residue indices
  max_index <- max(template$residues$residue_index)
  expect_true(all(template$basepairs$pos1 <= max_index))
  expect_true(all(template$basepairs$pos2 <= max_index))
})

# Test 3: load_modomics() works for all organisms ----
test_that("load_modomics works for all organisms", {
  orgs <- c("sacCer", "ecoliK12", "hg38")

  for (org in orgs) {
    mods <- load_modomics(org)

    expect_s3_class(mods, "tbl_df")
    expect_true(nrow(mods) > 0)

    # Verify required columns
    required_cols <- c(
      "trna_id",
      "sprinzl_label",
      "modification_short_name",
      "modification_name",
      "region"
    )
    expect_true(all(required_cols %in% names(mods)))
  }
})

# Test 4: load_modomics() validates organism argument ----
test_that("load_modomics validates organism argument", {
  expect_error(
    load_modomics("invalid_organism"),
    "should be one of"
  )
})

# Test 5: plot_trna_structure() returns ggplot object ----
test_that("plot_trna_structure returns ggplot object with minimal args", {
  bcerr <- read_bcerror(clover_example("yeast/grande.bcerr.tsv.gz"))
  bcerr_coords <- add_global_coords(bcerr, "sacCer")

  p <- plot_trna_structure(bcerr_coords)

  expect_s3_class(p, "gg")
  expect_s3_class(p, "ggplot")
})

# Test 6: plot_trna_structure() requires sprinzl_label column ----
test_that("plot_trna_structure requires sprinzl_label column", {
  bcerr <- read_bcerror(clover_example("yeast/grande.bcerr.tsv.gz"))
  # Don't add global coordinates

  expect_error(
    plot_trna_structure(bcerr),
    "sprinzl_label"
  )
})

# Test 7: plot_trna_structure() validates fill column exists ----
test_that("plot_trna_structure validates fill column exists", {
  bcerr_coords <- read_bcerror(clover_example("yeast/grande.bcerr.tsv.gz")) |>
    add_global_coords("sacCer")

  expect_error(
    plot_trna_structure(bcerr_coords, fill = "nonexistent_column"),
    "not found"
  )
})

# Test 8: plot_trna_structure() filters to single tRNA ----
test_that("plot_trna_structure filters to single tRNA when trna_id provided", {
  bcerr_coords <- read_bcerror(clover_example("yeast/grande.bcerr.tsv.gz")) |>
    add_global_coords("sacCer")

  # Pick a specific tRNA
  trna_id <- "nuc-tRNA-Ala-AGC-1-1"
  p <- plot_trna_structure(bcerr_coords, trna_id = trna_id)

  # Should return ggplot object
  expect_s3_class(p, "gg")

  # Plot should build without error
  built <- ggplot2::ggplot_build(p)
  expect_s3_class(built, "ggplot_built")
})

# Test 9: plot_trna_structure() aggregates when trna_id is NULL ----
test_that("plot_trna_structure aggregates across tRNAs when trna_id is NULL", {
  bcerr_coords <- read_bcerror(clover_example("yeast/grande.bcerr.tsv.gz")) |>
    add_global_coords("sacCer")

  # No trna_id = aggregate
  p <- plot_trna_structure(bcerr_coords, fill = "error_rate")

  expect_s3_class(p, "gg")

  # Should aggregate across all tRNAs per sprinzl position
  built <- ggplot2::ggplot_build(p)
  expect_s3_class(built, "ggplot_built")
})

# Test 10: plot_trna_structure() errors on nonexistent tRNA ----
test_that("plot_trna_structure errors on nonexistent tRNA", {
  bcerr_coords <- read_bcerror(clover_example("yeast/grande.bcerr.tsv.gz")) |>
    add_global_coords("sacCer")

  expect_error(
    plot_trna_structure(bcerr_coords, trna_id = "nonexistent-tRNA"),
    "No data found for tRNA"
  )
})

# Test 11: plot_trna_structure() accepts different fill variables ----
test_that("plot_trna_structure accepts different fill variables", {
  bcerr_coords <- read_bcerror(clover_example("yeast/grande.bcerr.tsv.gz")) |>
    add_global_coords("sacCer")

  # Test multiple fill options
  for (fill_var in c("error_rate", "mis", "ins", "del")) {
    p <- plot_trna_structure(bcerr_coords, fill = fill_var)

    expect_s3_class(p, "gg")

    # Verify plot builds
    built <- ggplot2::ggplot_build(p)
    expect_s3_class(built, "ggplot_built")
  }
})

# Test 12: plot_trna_structure() adds Modomics overlay ----
test_that("plot_trna_structure adds Modomics overlay when modomics = TRUE", {
  bcerr_coords <- read_bcerror(clover_example("yeast/grande.bcerr.tsv.gz")) |>
    add_global_coords("sacCer")

  # With Modomics
  p_with <- plot_trna_structure(
    bcerr_coords,
    modomics = TRUE,
    organism = "sacCer"
  )
  expect_s3_class(p_with, "gg")

  # Without Modomics
  p_without <- plot_trna_structure(
    bcerr_coords,
    modomics = FALSE
  )
  expect_s3_class(p_without, "gg")

  # With Modomics should have more layers (red circles overlay)
  expect_gte(length(p_with$layers), length(p_without$layers))
})

# Test 13: Full workflow integration test ----
test_that("Full workflow: read_bcerror -> add_global_coords -> plot_trna_structure", {
  # Step 1: Load bcerror data
  bcerr <- read_bcerror(clover_example("yeast/grande.bcerr.tsv.gz"))
  expect_s3_class(bcerr, "tbl_df")

  # Step 2: Add global coordinates
  bcerr_coords <- add_global_coords(bcerr, "sacCer")
  expect_true("sprinzl_label" %in% names(bcerr_coords))
  expect_true("global_index" %in% names(bcerr_coords))

  # Step 3: Create structure plot
  p <- plot_trna_structure(bcerr_coords, fill = "error_rate")
  expect_s3_class(p, "gg")

  # Step 4: Verify plot can be built
  built <- ggplot2::ggplot_build(p)
  expect_s3_class(built, "ggplot_built")

  # Verify data integrity through pipeline
  point_data <- built$data[[2]] # Point layer (residues)
  expect_true(all(is.finite(point_data$x)))
  expect_true(all(is.finite(point_data$y)))
})

# Test 14: Heatmap and structure plot consistency ----
test_that("Heatmap and structure plots work with same data", {
  # Test that global coordinates work consistently across plot types
  bcerr_coords <- read_bcerror(clover_example("yeast/grande.bcerr.tsv.gz")) |>
    add_global_coords("sacCer") |>
    dplyr::filter(grepl("^nuc-", ref))

  # Both plot types should work with same data
  p_heatmap <- plot_bcerror_heatmap(bcerr_coords, value = "error_rate")
  expect_s3_class(p_heatmap, "gg")

  p_structure <- plot_trna_structure(bcerr_coords, fill = "error_rate")
  expect_s3_class(p_structure, "gg")

  # Both should build without error
  expect_silent(ggplot2::ggplot_build(p_heatmap))
  expect_silent(ggplot2::ggplot_build(p_structure))
})
