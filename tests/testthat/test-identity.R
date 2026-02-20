test_that("identity_organisms returns all 3 organisms", {
  orgs <- identity_organisms()
  expect_type(orgs, "character")
  expect_length(orgs, 3)
  expect_contains(orgs, "Escherichia coli")
  expect_contains(orgs, "Saccharomyces cerevisiae")
  expect_contains(orgs, "Homo sapiens")
})

test_that("identity_elements returns expected columns and types", {
  result <- identity_elements("Escherichia coli")

  expect_s3_class(result, "tbl_df")
  expect_named(
    result,
    c(
      "amino_acid",
      "domain",
      "sprinzl_pos",
      "nucleotide",
      "region",
      "type",
      "strength",
      "pair_pos",
      "pair_type",
      "against_aars",
      "universal",
      "description"
    )
  )

  expect_type(result$amino_acid, "character")
  expect_type(result$domain, "character")
  expect_type(result$sprinzl_pos, "integer")
  expect_type(result$nucleotide, "character")
  expect_type(result$region, "character")
  expect_type(result$type, "character")
  expect_type(result$strength, "character")
  expect_type(result$pair_pos, "integer")
  expect_type(result$pair_type, "character")
  expect_type(result$against_aars, "character")
  expect_type(result$universal, "logical")
  expect_type(result$description, "character")
})

test_that("identity_elements filters by amino_acid", {
  result <- identity_elements("Escherichia coli", amino_acid = "Ala")
  expect_equal(unique(result$amino_acid), "Ala")
  expect_gt(nrow(result), 0)
})

test_that("identity_elements filters by type", {
  dets <- identity_elements("Escherichia coli", type = "determinant")
  expect_equal(unique(dets$type), "determinant")

  antis <- identity_elements("Escherichia coli", type = "antideterminant")
  expect_equal(unique(antis$type), "antideterminant")

  both <- identity_elements("Escherichia coli", type = "both")
  expect_equal(sort(unique(both$type)), c("antideterminant", "determinant"))
})

test_that("Ala G3-U70 wobble pair is present for Bacteria", {
  result <- identity_elements(
    "Escherichia coli",
    amino_acid = "Ala",
    type = "determinant"
  )

  g3 <- result[result$sprinzl_pos == 3 & !is.na(result$sprinzl_pos), ]
  expect_equal(nrow(g3), 1)
  expect_equal(g3$nucleotide, "G")
  expect_equal(g3$pair_pos, 70L)
  expect_equal(g3$pair_type, "wobble")
  expect_true(g3$universal)
})

test_that("identity_elements errors for unsupported organism", {
  expect_snapshot(
    identity_elements("Mus musculus"),
    error = TRUE
  )
})

test_that("Eukarya organisms return identical elements", {
  sc <- identity_elements("Saccharomyces cerevisiae")
  hs <- identity_elements("Homo sapiens")
  expect_equal(sc, hs)
})

test_that("all 20 amino acids have Bacteria determinants", {
  result <- identity_elements(
    "Escherichia coli",
    type = "determinant"
  )
  aa_list <- sort(unique(result$amino_acid))
  expect_length(aa_list, 20)
})

test_that("map_identity_to_trna converts Sprinzl to seq positions", {
  elements <- identity_elements(
    "Escherichia coli",
    amino_acid = "Ala",
    type = "determinant"
  )

  # Create mock Sprinzl coordinates
  mock_coords <- dplyr::tibble(
    trna_id = "tRNA-Ala-GGC",
    pos = c(1L, 2L, 3L, 68L, 69L, 70L, 71L, 73L),
    sprinzl_label = c("1", "2", "3", "68", "69", "70", "71", "73"),
    global_index = seq_len(8),
    region = "acceptor_stem",
    residue = "N"
  )

  result <- map_identity_to_trna(elements, mock_coords, "tRNA-Ala-GGC")

  expect_s3_class(result, "tbl_df")
  expect_true("pos" %in% names(result))

  # G3-U70 should map; G-1 and G20 should be dropped
  mapped_sprinzl <- result$sprinzl_pos
  expect_true(3L %in% mapped_sprinzl)
  expect_true(70L %in% mapped_sprinzl)
  expect_false(-1L %in% mapped_sprinzl)
})

test_that("map_identity_to_trna errors for missing tRNA", {
  elements <- identity_elements("Escherichia coli", amino_acid = "Ala")
  mock_coords <- dplyr::tibble(
    trna_id = "other-tRNA",
    pos = 1L,
    sprinzl_label = "1",
    global_index = 1L,
    region = "acceptor_stem",
    residue = "G"
  )

  expect_snapshot(
    map_identity_to_trna(elements, mock_coords, "tRNA-Ala-GGC"),
    error = TRUE
  )
})
