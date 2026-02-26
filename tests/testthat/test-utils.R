test_that("shorten_trna_names strips prefix, tRNA-, and copy suffix", {
  expect_equal(
    shorten_trna_names("host-tRNA-Ser-CGT-1-1"),
    "Ser-CGT"
  )
  expect_equal(
    shorten_trna_names("phage-tRNA-Glu-TTC-2-1"),
    "Glu-TTC"
  )
})

test_that("shorten_trna_names works without source prefix", {
  expect_equal(shorten_trna_names("tRNA-Ala-AGC-3-1"), "Ala-AGC")
})

test_that("shorten_trna_names handles vector input", {
  input <- c(
    "host-tRNA-Ser-CGT-1-1",
    "phage-tRNA-Glu-TTC-2-1",
    "tRNA-Ala-AGC-3-1"
  )
  expect_equal(
    shorten_trna_names(input),
    c("Ser-CGT", "Glu-TTC", "Ala-AGC")
  )
})

test_that("shorten_trna_names respects custom strip_prefix", {
  expect_equal(
    shorten_trna_names("custom-tRNA-Ser-CGT-1-1", strip_prefix = "^custom-"),
    "Ser-CGT"
  )
  # Default pattern should not strip "custom-"
  expect_equal(
    shorten_trna_names("custom-tRNA-Ser-CGT-1-1"),
    "custom-tRNA-Ser-CGT"
  )
})

test_that("shorten_trna_names passes through non-matching strings", {
  expect_equal(shorten_trna_names("some-other-name"), "some-other-name")
})

test_that("dna_to_rna_anticodon converts T to U in anticodon", {
  expect_equal(
    dna_to_rna_anticodon("tRNA-Glu-TTC-1-1"),
    "tRNA-Glu-UUC-1-1"
  )
})

test_that("dna_to_rna_anticodon passes through non-matching strings", {
  expect_equal(dna_to_rna_anticodon("other"), "other")
})
