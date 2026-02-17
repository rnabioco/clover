# Mock modification dictionary for unit tests
mock_mod_dict <- function() {
  dplyr::tibble(
    new_abbrev = c("D", "P", "K", "t"),
    short_name = c("D", "pseudoU", "m1G", "m5U"),
    name = c(
      "dihydrouridine",
      "pseudouridine",
      "1-methylguanosine",
      "5-methyluridine"
    ),
    reference_moiety = c("U", "U", "G", "U")
  )
}

test_that(".extract_mod_positions finds known modifications", {
  mod_dict <- mock_mod_dict()
  # Positions:    1234567
  seq <- "AUGDCGP"

  result <- .extract_mod_positions(seq, mod_dict)

  expect_s3_class(result, "tbl_df")
  expect_named(result, c("pos", "mod_full", "mod1"))
  expect_equal(result$pos, c(4L, 7L))
  expect_equal(result$mod1, c("D", "pseudoU"))
  expect_equal(
    result$mod_full,
    c("dihydrouridine", "pseudouridine")
  )
})

test_that(".extract_mod_positions returns empty for unmodified sequence", {
  mod_dict <- mock_mod_dict()
  result <- .extract_mod_positions("AUGCAUGC", mod_dict)

  expect_equal(nrow(result), 0)
  expect_named(result, c("pos", "mod_full", "mod1"))
})

test_that(".extract_mod_positions handles empty/NA input", {
  mod_dict <- mock_mod_dict()
  expect_equal(nrow(.extract_mod_positions("", mod_dict)), 0)
  expect_equal(nrow(.extract_mod_positions(NA, mod_dict)), 0)
})

test_that(".strip_modifications replaces mod chars with parent base", {
  mod_dict <- mock_mod_dict()
  # D -> U, P -> U
  result <- .strip_modifications("AUGDCGP", mod_dict)
  expect_equal(result, "AUGUCGU")
})

test_that(".strip_modifications replaces unknowns with N", {
  mod_dict <- mock_mod_dict()
  # Z is not in mod_dict and not a standard base
  result <- .strip_modifications("AUGZCG", mod_dict)
  expect_equal(result, "AUGNCG")
})

test_that(".strip_modifications handles empty/NA input", {
  mod_dict <- mock_mod_dict()
  expect_equal(.strip_modifications("", mod_dict), "")
  expect_equal(.strip_modifications(NA, mod_dict), "")
})

test_that(".transfer_positions maps through alignment", {
  skip_if_not_installed("pwalign")
  pattern <- Biostrings::DNAString("ACGTACGT")
  subject <- Biostrings::DNAString("ACGTACGT")
  aln <- pwalign::pairwiseAlignment(
    pattern,
    subject,
    type = "local"
  )

  mod_positions <- dplyr::tibble(
    pos = c(2L, 5L),
    mod_full = c("mod_a", "mod_b"),
    mod1 = c("mA", "mB")
  )

  result <- .transfer_positions(aln, mod_positions)

  expect_equal(nrow(result), 2)
  expect_equal(result$pos, c(2L, 5L))
  expect_equal(result$mod1, c("mA", "mB"))
})

test_that(".transfer_positions handles gaps in alignment", {
  skip_if_not_installed("pwalign")
  # Pattern has an insertion relative to subject (use global
  # alignment so all positions are covered)
  pattern <- Biostrings::DNAString("ACGTTACGT")
  subject <- Biostrings::DNAString("ACGTACGT")
  aln <- pwalign::pairwiseAlignment(
    pattern,
    subject,
    type = "global"
  )

  mod_positions <- dplyr::tibble(
    pos = c(1L, 9L),
    mod_full = c("mod_a", "mod_b"),
    mod1 = c("mA", "mB")
  )

  result <- .transfer_positions(aln, mod_positions)

  # Position 1 should map to 1, position 9 to 8
  expect_equal(result$pos[result$mod1 == "mA"], 1L)
  expect_equal(result$pos[result$mod1 == "mB"], 8L)
})

test_that(".find_aa_candidates matches amino acid in FASTA names", {
  fasta_names <- c(
    "tRNA-Ala-AGC-1",
    "tRNA-Ala-TGC-2",
    "tRNA-Gly-GCC-1",
    "tRNA-Met-CAT-1"
  )

  ala_idx <- .find_aa_candidates("Ala", fasta_names)
  expect_equal(sort(ala_idx), c(1L, 2L))

  gly_idx <- .find_aa_candidates("Gly", fasta_names)
  expect_equal(gly_idx, 3L)
})

test_that(".find_aa_candidates maps Ini to Met variants", {
  fasta_names <- c(
    "tRNA-Met-CAT-1",
    "tRNA-iMet-CAT-1",
    "tRNA-Ala-AGC-1"
  )

  ini_idx <- .find_aa_candidates("Ini", fasta_names)
  expect_true(1L %in% ini_idx)
  expect_true(2L %in% ini_idx)
  expect_false(3L %in% ini_idx)
})

test_that(".find_aa_candidates returns empty for unknown AA", {
  fasta_names <- c("tRNA-Ala-AGC-1", "tRNA-Gly-GCC-1")
  result <- .find_aa_candidates("Xyz", fasta_names)
  expect_length(result, 0)
})

test_that("fetch_modomics_mods returns expected structure with mocked API", {
  skip_if_not_installed("httr2")
  skip_if_not_installed("jsonlite")
  skip_if_not_installed("pwalign")

  mock_mods_json <- jsonlite::toJSON(
    list(
      "1" = list(
        new_abbrev = "D",
        short_name = "D",
        name = "dihydrouridine",
        reference_moiety = list("U")
      ),
      "2" = list(
        new_abbrev = "P",
        short_name = "pseudoU",
        name = "pseudouridine",
        reference_moiety = list("U")
      )
    ),
    auto_unbox = TRUE
  )

  mock_seqs_json <- jsonlite::toJSON(
    list(
      "1" = list(
        subtype = "Ala",
        anticodon = "AGC",
        seq = "AUGDCGAUGPCGA"
      )
    ),
    auto_unbox = TRUE
  )

  local_mocked_bindings(
    .fetch_modomics_modifications = function(cache_dir = NULL) {
      entries <- jsonlite::fromJSON(
        as.character(mock_mods_json),
        simplifyVector = FALSE
      )
      dplyr::tibble(
        new_abbrev = vapply(
          entries,
          \(e) e$new_abbrev %||% NA_character_,
          character(1)
        ),
        short_name = vapply(
          entries,
          \(e) e$short_name %||% NA_character_,
          character(1)
        ),
        name = vapply(
          entries,
          \(e) e$name %||% NA_character_,
          character(1)
        ),
        reference_moiety = vapply(
          entries,
          function(e) {
            rm_val <- e$reference_moiety
            if (is.null(rm_val) || length(rm_val) == 0) {
              return(NA_character_)
            }
            rm_val[[1]]
          },
          character(1)
        )
      )
    },
    .fetch_modomics_sequences = function(organism, cache_dir = NULL) {
      entries <- jsonlite::fromJSON(
        as.character(mock_seqs_json),
        simplifyVector = FALSE
      )
      dplyr::tibble(
        subtype = vapply(
          entries,
          \(e) e$subtype %||% NA_character_,
          character(1)
        ),
        anticodon = vapply(
          entries,
          \(e) e$anticodon %||% NA_character_,
          character(1)
        ),
        seq = vapply(
          entries,
          \(e) e$seq %||% NA_character_,
          character(1)
        )
      )
    }
  )

  # Create a minimal FASTA with a matching sequence
  # AUGDCGAUGPCGA stripped = AUGUCGAUGUCGA, DNA = ATGTCGATGTCGA
  ref_seq <- Biostrings::DNAStringSet("ATGTCGATGTCGA")
  names(ref_seq) <- "tRNA-Ala-AGC-1"

  result <- fetch_modomics_mods(ref_seq, "test organism")

  expect_s3_class(result, "tbl_df")
  expect_named(result, c("ref", "pos", "mod_full", "mod1"))
  expect_gt(nrow(result), 0)
  expect_equal(result$ref[1], "tRNA-Ala-AGC-1")
})

test_that("fetch_modomics_mods works with yeast data", {
  skip_if_not_installed("httr2")
  skip_if_not_installed("jsonlite")
  skip_if_not_installed("pwalign")
  skip_if_offline()

  fa <- clover_example("yeast/trna-ref.fa.gz")
  mods <- fetch_modomics_mods(fa, "Saccharomyces cerevisiae")

  expect_s3_class(mods, "tbl_df")
  expect_named(mods, c("ref", "pos", "mod_full", "mod1"))
  expect_gt(nrow(mods), 0)
})
