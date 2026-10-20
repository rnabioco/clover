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

test_that("normalize_trna_id maps every source spelling to one form", {
  x <- c(
    "host-tRNA-Ala-GGC-1-1",
    "phage-tRNA-Arg-TCT",
    "nuc-tRNA-Ala-AGC-1-1",
    "mito-tRNA-fMet-CAU",
    "tRNA-Ile-GAU-1-1",
    "trna-Gln-TTG",
    "tRNA-Gly-GCC-1-1-uncharged",
    "tRNA-Ile2-CAT-1-1",
    "tRNA-Ala-GGC"
  )
  expect_equal(
    normalize_trna_id(x),
    c(
      "tRNA-Ala-GGC",
      "tRNA-Arg-TCT",
      "tRNA-Ala-AGC",
      "tRNA-fMet-CAT",
      "tRNA-Ile-GAT",
      "tRNA-Gln-TTG",
      "tRNA-Gly-GCC",
      "tRNA-Ile2-CAT",
      "tRNA-Ala-GGC"
    )
  )
})

test_that("normalize_trna_id keeps the copy suffix only when asked", {
  expect_equal(
    normalize_trna_id("phage-tRNA-Thr-CGT-2-1", copy = TRUE),
    "tRNA-Thr-CGT-2-1"
  )
  expect_equal(
    normalize_trna_id("tRNA-Gly-GCC-1-1-uncharged", copy = TRUE),
    "tRNA-Gly-GCC-1-1"
  )
  expect_equal(
    normalize_trna_id("mito-tRNA-fMet-CAU", copy = TRUE),
    "tRNA-fMet-CAT"
  )
})

test_that("normalize_trna_id is idempotent and returns NA for non-tRNA IDs", {
  x <- c("host-tRNA-Ala-GGC-1-1", "tRNA-Ile-GAU-1-1", "mito-tRNA-fMet-CAU")
  once <- normalize_trna_id(x)
  expect_equal(normalize_trna_id(once), once)
  expect_equal(
    normalize_trna_id(normalize_trna_id(x, copy = TRUE), copy = TRUE),
    normalize_trna_id(x, copy = TRUE)
  )
  expect_equal(
    normalize_trna_id(c("tRNA-Ala", "not a tRNA", NA)),
    rep(NA_character_, 3)
  )
  expect_equal(normalize_trna_id(character()), character())
})

test_that("every bundled Sprinzl coordinate ID resolves", {
  files <- list.files(
    clover_example("sprinzl"),
    pattern = "coords\\.tsv\\.gz$",
    full.names = TRUE
  )
  expect_gt(length(files), 0)
  for (f in files) {
    ids <- unique(read_sprinzl_coords(f)$trna_id)
    expect_equal(
      ids[is.na(normalize_trna_id(ids))],
      character(),
      info = basename(f)
    )
  }
})

test_that("bundled E. coli FASTA, coordinates and bcerror share keys", {
  fa_names <- names(Biostrings::readDNAStringSet(
    clover_example("ecoli/trna_only.fa.gz")
  ))
  coords <- read_sprinzl_coords(
    clover_example("sprinzl/ecoliK12_global_coords.tsv.gz")
  )
  bc <- read_bcerror(
    clover_example(
      "ecoli/summary/tables/wt-15-ctl-01/wt-15-ctl-01.bcerror.tsv.gz"
    )
  )

  expect_equal(fa_names[is.na(normalize_trna_id(fa_names))], character())
  expect_equal(
    unique(bc$ref)[is.na(normalize_trna_id(unique(bc$ref)))],
    character()
  )

  coord_keys <- unique(normalize_trna_id(coords$trna_id))
  expect_gt(length(intersect(normalize_trna_id(bc$ref), coord_keys)), 0)
  expect_gt(length(intersect(normalize_trna_id(fa_names), coord_keys)), 0)
  # Same key, with and without the charging-reference variant
  expect_equal(
    normalize_trna_id("host-tRNA-Gly-GCC-1-1-uncharged"),
    normalize_trna_id("tRNA-Gly-GCC-1-1")
  )
})
