test_that("aa_colors returns named character vector", {
  colors <- aa_colors()
  expect_type(colors, "character")
  expect_true(all(nzchar(names(colors))))
  expect_length(colors, 23)
})

test_that("aa_colors includes standard amino acids", {
  colors <- aa_colors()
  standard <- c(
    "Ala",
    "Arg",
    "Asn",
    "Asp",
    "Cys",
    "Gln",
    "Glu",
    "Gly",
    "His",
    "Ile",
    "Leu",
    "Lys",
    "Met",
    "Phe",
    "Pro",
    "Ser",
    "Thr",
    "Trp",
    "Tyr",
    "Val"
  )
  expect_true(all(standard %in% names(colors)))
})

test_that("aa_colors includes special amino acids", {
  colors <- aa_colors()
  expect_true(all(c("fMet", "Ile2", "SeC") %in% names(colors)))
})

test_that("charging_colors returns named character vector", {
  colors <- charging_colors()
  expect_type(colors, "character")
  expect_named(colors, c("Charged", "Uncharged"))
  expect_length(colors, 2)
})
