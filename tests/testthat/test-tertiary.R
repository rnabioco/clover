test_that("tertiary_contacts returns expected columns and types", {
  result <- tertiary_contacts("Escherichia coli")

  expect_s3_class(result, "tbl_df")
  expect_named(
    result,
    c(
      "contact_id",
      "contact_type",
      "domain",
      "pos1",
      "pos2",
      "pos3",
      "nuc1",
      "nuc2",
      "nuc3",
      "interaction",
      "universal",
      "description"
    )
  )

  expect_type(result$contact_id, "character")
  expect_type(result$contact_type, "character")
  expect_type(result$pos1, "integer")
  expect_type(result$pos2, "integer")
  expect_type(result$pos3, "integer")
  expect_type(result$universal, "logical")
})

test_that("tertiary_contacts returns 9 rows for Bacteria", {
  result <- tertiary_contacts("Escherichia coli")
  expect_equal(nrow(result), 9)
  expect_equal(unique(result$domain), "Bacteria")
})

test_that("tertiary_contacts splits 4 triples and 5 pairs", {
  result <- tertiary_contacts("Escherichia coli")
  expect_equal(sum(result$contact_type == "triple"), 4)
  expect_equal(sum(result$contact_type == "pair"), 5)

  pairs <- result[result$contact_type == "pair", ]
  expect_all_true(is.na(pairs$pos3))
  expect_all_true(is.na(pairs$nuc3))
})

test_that("U8-A14-A21 has expected positions and bases", {
  result <- tertiary_contacts("Escherichia coli")
  triple <- result[result$contact_id == "U8-A14-A21", ]

  expect_equal(nrow(triple), 1)
  expect_equal(triple$contact_type, "triple")
  expect_equal(triple$pos1, 8L)
  expect_equal(triple$pos2, 14L)
  expect_equal(triple$pos3, 21L)
  expect_equal(triple$nuc1, "U")
  expect_equal(triple$nuc2, "A")
  expect_equal(triple$nuc3, "A")
  expect_true(triple$universal)
})

test_that("G15-C48 Levitt pair is encoded as a pair", {
  result <- tertiary_contacts("Escherichia coli")
  levitt <- result[result$contact_id == "G15-C48", ]

  expect_equal(nrow(levitt), 1)
  expect_equal(levitt$contact_type, "pair")
  expect_equal(levitt$pos1, 15L)
  expect_equal(levitt$pos2, 48L)
  expect_equal(levitt$nuc1, "G")
  expect_equal(levitt$nuc2, "C")
  expect_true(is.na(levitt$pos3))
  expect_equal(levitt$interaction, "Levitt")
})

test_that("tertiary_contacts works for all supported organisms", {
  for (org in identity_organisms()) {
    result <- tertiary_contacts(org)
    expect_equal(nrow(result), 9)
  }
})

test_that("tertiary_contacts errors for unsupported organism", {
  expect_snapshot(
    tertiary_contacts("Mus musculus"),
    error = TRUE
  )
})
