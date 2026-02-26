test_that("read_fasta returns DNAStringSet with correct names", {
  fa <- clover_example("ecoli/trna_only.fa.gz")
  result <- read_fasta(fa)

  expect_s4_class(result, "DNAStringSet")
  expect_true(length(result) > 0)
  expect_true(all(nchar(as.character(result)) > 0))
})

test_that("read_mod_annotations returns tibble with expected columns", {
  tmp <- withr::local_tempfile(fileext = ".tsv")
  writeLines(
    c("ref\tpos\tmod_full\tmod1", "tRNA-Ala\t20\t1-methyladenosine\tm1A"),
    tmp
  )
  result <- read_mod_annotations(tmp)

  expect_s3_class(result, "tbl_df")
  expect_named(result, c("ref", "pos", "mod_full", "mod1"))
  expect_equal(nrow(result), 1)
})

test_that("create_clover builds SE from ecoli test data", {
  config_path <- clover_example("ecoli/config.yaml")

  se <- create_clover(config_path, types = c("charging", "bcerror"))

  expect_s4_class(se, "SummarizedExperiment")

  # Check counts assay exists and has correct dimensions
  counts <- SummarizedExperiment::assay(se, "counts")
  expect_true(is.matrix(counts))
  expect_equal(ncol(counts), 6) # 6 samples
  expect_true(nrow(counts) > 0)

  # Column names should be sample IDs
  expect_true(all(grepl("^wt-15-", colnames(counts))))

  # Check colData
  cd <- as.data.frame(SummarizedExperiment::colData(se))
  expect_true("sample_id" %in% names(cd))
  expect_equal(nrow(cd), 6)

  # Check metadata has charging and bcerror
  meta <- S4Vectors::metadata(se)
  expect_true("config" %in% names(meta))
  expect_s3_class(meta$charging, "tbl_df")
  expect_true(nrow(meta$charging) > 0)
  expect_s3_class(meta$bcerror, "tbl_df")
  expect_true(nrow(meta$bcerror) > 0)

  # Check rowData
  rd <- as.data.frame(SummarizedExperiment::rowData(se))
  expect_true("ref" %in% names(rd))
})

test_that("create_clover with sample_info merges metadata", {
  config_path <- clover_example("ecoli/config.yaml")

  sample_info <- data.frame(
    sample_id = c(
      "wt-15-ctl-01",
      "wt-15-ctl-02",
      "wt-15-ctl-03",
      "wt-15-inf-01",
      "wt-15-inf-02",
      "wt-15-inf-03"
    ),
    condition = rep(c("ctl", "inf"), each = 3)
  )

  se <- create_clover(
    config_path,
    types = "charging",
    sample_info = sample_info
  )

  cd <- as.data.frame(SummarizedExperiment::colData(se))
  expect_true("condition" %in% names(cd))
  expect_equal(sort(unique(cd$condition)), c("ctl", "inf"))
})

test_that("create_clover loads FASTA into metadata", {
  config_path <- clover_example("ecoli/config.yaml")

  se <- create_clover(config_path, types = "charging")

  meta <- S4Vectors::metadata(se)
  expect_true("fasta" %in% names(meta))
  expect_s4_class(meta$fasta, "DNAStringSet")
})
