# Test coordinate loading and auto-classification

test_that("available_organisms returns expected organisms", {
  orgs <- available_organisms()

  expect_true("ecoliK12" %in% orgs)
  expect_true("sacCer" %in% orgs)
  expect_true("hg38" %in% orgs)
})

test_that("available_coord_groups lists groups for organism", {
  groups <- available_coord_groups("sacCer")

  expect_s3_class(groups, "tbl_df")
  expect_true(all(c("organism", "offset", "type", "n_trnas", "file") %in% names(groups)))
  expect_true(nrow(groups) >= 1)

  # Should have type1 and type2 groups
  expect_true("type1" %in% groups$type)
  expect_true("type2" %in% groups$type)
})

test_that("load_global_coords loads all groups by default", {
  coords <- load_global_coords("sacCer")

  expect_s3_class(coords, "tbl_df")
  expect_true(nrow(coords) > 0)
  expect_true(all(c("trna_id", "seq_index", "global_index", "region", "offset", "type") %in% names(coords)))

  # Should have multiple groups
  groups <- unique(paste0(coords$offset, "_", coords$type))
  expect_true(length(groups) > 1)

  # Check that global_index is positive integers
  expect_true(all(coords$global_index > 0))

  # Check that regions are present
  expect_true("D-loop" %in% coords$region)
  expect_true("anticodon-loop" %in% coords$region)
})

test_that("load_global_coords can filter by type", {
  coords_type1 <- load_global_coords("sacCer", type = "type1")
  coords_type2 <- load_global_coords("sacCer", type = "type2")

  expect_true(all(coords_type1$type == "type1"))
  expect_true(all(coords_type2$type == "type2"))

  # Type II should have Leu, Ser, Tyr tRNAs
  type2_trnas <- unique(coords_type2$trna_id)
  expect_true(any(grepl("Leu|Ser|Tyr", type2_trnas)))
})

test_that("load_global_coords can filter by type and offset", {
  coords <- load_global_coords("sacCer", type = "type1", offset = 0)

  expect_true(all(coords$type == "type1"))
  expect_true(all(coords$offset == 0))
})

test_that("add_global_coords auto-classifies tRNAs", {
  bcerr_path <- clover_example("yeast/grande.bcerr.tsv.gz")
  bcerr <- read_bcerror(bcerr_path)

  # Expect warning about unmatched tRNAs (mito, etc)
  expect_warning(
    bcerr_coords <- add_global_coords(bcerr, "sacCer"),
    "not found in coordinate groups"
  )

  expect_true("global_index" %in% names(bcerr_coords))
  expect_true("sprinzl_label" %in% names(bcerr_coords))
  expect_true("region" %in% names(bcerr_coords))
  expect_true("offset" %in% names(bcerr_coords))
  expect_true("type" %in% names(bcerr_coords))

  # Nuclear tRNAs should have coordinates
  nuc_bcerr <- dplyr::filter(bcerr_coords, grepl("^nuc-", ref))
  expect_true(sum(!is.na(nuc_bcerr$global_index)) > 0)
})

test_that("sprinzl position 55 aligns correctly within groups", {
  bcerr_path <- clover_example("yeast/grande.bcerr.tsv.gz")
  bcerr <- read_bcerror(bcerr_path)

  bcerr_coords <- suppressWarnings(add_global_coords(bcerr, "sacCer"))

  # Within each group, position 55 should map to exactly one global_index
  pos55 <- bcerr_coords |>
    dplyr::filter(sprinzl_label == "55", !is.na(global_index)) |>
    dplyr::group_by(type, offset) |>
    dplyr::summarize(n_indices = dplyr::n_distinct(global_index), .groups = "drop")

  expect_true(all(pos55$n_indices == 1),
              label = "Position 55 should map to single global_index per group")
})

test_that("get_global_labels returns named vector for single group", {
  coords <- load_global_coords("sacCer", type = "type1", offset = 0)
  labels <- get_global_labels(coords)

  expect_type(labels, "character")
  expect_true(length(labels) > 0)
  expect_true(!is.null(names(labels)))

  # Should include common positions
  expect_true("55" %in% labels)
  expect_true("34" %in% labels)
})

test_that("get_global_labels returns list for multiple groups", {
  coords <- load_global_coords("sacCer")
  labels <- get_global_labels(coords)

  # Should be a list when multiple groups
  expect_type(labels, "list")
  expect_true(length(labels) > 1)
})

test_that("get_region_bounds returns region boundaries", {
  coords <- load_global_coords("sacCer", type = "type1", offset = 0)
  regions <- get_region_bounds(coords)

  expect_s3_class(regions, "tbl_df")
  expect_true(all(c("region", "start", "end") %in% names(regions)))
  expect_true(all(regions$start <= regions$end))
})

test_that("adapter offset correctly aligns bcerror with coordinates", {
  bcerr <- read_bcerror(clover_example("yeast/grande.bcerr.tsv.gz"))
  bcerr_coords <- suppressWarnings(add_global_coords(bcerr, "sacCer"))

  # Check that is_adapter column was added
  expect_true("is_adapter" %in% names(bcerr_coords))

  # Nuclear tRNA positions 25-97 should have coordinate matches
  nuc_tRNA <- bcerr_coords |>
    dplyr::filter(grepl("^nuc-", ref), pos >= 25, pos <= 97)

  matched <- sum(!is.na(nuc_tRNA$global_index))
  total <- nrow(nuc_tRNA)
  match_rate <- matched / total

  expect_gt(match_rate, 0.80,
            label = paste("Match rate:", round(match_rate * 100, 1), "%"))

  # Adapter positions (1-24) should have NO matches
  adapter_region <- bcerr_coords |>
    dplyr::filter(grepl("^nuc-", ref), pos <= 24)

  expect_true(all(is.na(adapter_region$global_index)),
              label = "Adapter positions should not match coordinates")

  # is_adapter flag should be TRUE for adapter positions
  expect_true(all(adapter_region$is_adapter),
              label = "is_adapter should be TRUE for positions 1-24")

  # Most tail positions (beyond standard tRNA length) should have NO matches
  # Note: Type II tRNAs (Leu, Ser, Tyr) may have valid coords up to ~97+24=121
  far_tail_region <- bcerr_coords |>
    dplyr::filter(grepl("^nuc-", ref), pos >= 110)

  # Far tail should have very few or no matches
  far_tail_match_rate <- sum(!is.na(far_tail_region$global_index)) / nrow(far_tail_region)
  expect_lt(far_tail_match_rate, 0.05,
            label = "Far tail positions (110+) should have <5% match rate")
})

test_that("classify_trna_type identifies Type I and Type II", {
  type1 <- c("tRNA-Ala-GGC-1-1", "nuc-tRNA-Arg-ACG-1-1", "tRNA-Gly-CCC-1")
  type2 <- c("tRNA-Leu-CAA-1-1", "nuc-tRNA-Ser-GCT-1-1", "tRNA-Tyr-GTA-1")

  expect_equal(classify_trna_type(type1), rep("Type I", 3))
  expect_equal(classify_trna_type(type2), rep("Type II", 3))
})

test_that("mitochondrial tRNAs are not in coordinate groups", {
  bcerr <- read_bcerror(clover_example("yeast/grande.bcerr.tsv.gz"))
  bcerr_coords <- suppressWarnings(add_global_coords(bcerr, "sacCer"))

  # Mitochondrial tRNAs should have NA coordinates

  mito <- bcerr_coords |>
    dplyr::filter(grepl("^mito-", ref))

  if (nrow(mito) > 0) {
    expect_true(all(is.na(mito$global_index)),
                label = "Mitochondrial tRNAs should not have coordinates")
  }
})
