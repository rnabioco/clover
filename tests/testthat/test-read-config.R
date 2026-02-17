test_that("read_pipeline_config parses YAML with inline samples", {
  tmp <- normalizePath(withr::local_tempdir())

  config_text <- paste(
    "samples:",
    "  sample1: data/sample1.pod5",
    "  sample2: data/sample2.pod5",
    "output_dir: results",
    "fasta: ref/trna.fa",
    sep = "\n"
  )
  config_path <- file.path(tmp, "config.yaml")
  writeLines(config_text, config_path)

  cfg <- read_pipeline_config(config_path)

  expect_type(cfg, "list")
  expect_named(cfg, c("samples", "output_dir", "fasta", "config_dir"))
  expect_s3_class(cfg$samples, "tbl_df")
  expect_equal(nrow(cfg$samples), 2)
  expect_named(cfg$samples, c("sample_id", "data_path"))
  expect_equal(cfg$samples$sample_id, c("sample1", "sample2"))

  # Paths should be resolved relative to config dir
  expect_true(startsWith(cfg$output_dir, tmp))
  expect_true(startsWith(cfg$fasta, tmp))
})

test_that("read_pipeline_config parses YAML with TSV sample file (headered)", {
  tmp <- withr::local_tempdir()

  # Create a samples TSV with header
  samples_tsv <- file.path(tmp, "samples.tsv")
  writeLines(
    c("sample_id\tdata_path", "s1\tdata/s1.pod5", "s2\tdata/s2.pod5"),
    samples_tsv
  )

  config_text <- paste(
    "samples: samples.tsv",
    "output_dir: results",
    "fasta: ref/trna.fa",
    sep = "\n"
  )
  config_path <- file.path(tmp, "config.yaml")
  writeLines(config_text, config_path)

  cfg <- read_pipeline_config(config_path)

  expect_s3_class(cfg$samples, "tbl_df")
  expect_equal(nrow(cfg$samples), 2)
  expect_true("sample_id" %in% names(cfg$samples))
})

test_that("read_pipeline_config handles headerless samples.tsv", {
  tmp <- withr::local_tempdir()

  # Create a headerless samples TSV (real pipeline format)
  samples_tsv <- file.path(tmp, "samples.tsv")
  writeLines(
    c("s1\t/data/s1.pod5", "s2\t/data/s2.pod5"),
    samples_tsv
  )

  config_text <- paste(
    "samples: samples.tsv",
    "output_dir: results",
    "fasta: ref/trna.fa",
    sep = "\n"
  )
  config_path <- file.path(tmp, "config.yaml")
  writeLines(config_text, config_path)

  cfg <- read_pipeline_config(config_path)

  expect_s3_class(cfg$samples, "tbl_df")
  expect_equal(nrow(cfg$samples), 2)
  expect_named(cfg$samples, c("sample_id", "data_path"))
  expect_equal(cfg$samples$sample_id, c("s1", "s2"))
})

test_that("read_pipeline_config supports output_directory key", {
  tmp <- normalizePath(withr::local_tempdir())

  config_text <- paste(
    "samples:",
    "  sample1: data/s1.pod5",
    "output_directory: results",
    "fasta: ref/trna.fa",
    sep = "\n"
  )
  config_path <- file.path(tmp, "config.yaml")
  writeLines(config_text, config_path)

  cfg <- read_pipeline_config(config_path)

  expect_true(startsWith(cfg$output_dir, tmp))
  expect_true(grepl("results$", cfg$output_dir))
})

test_that("read_pipeline_config resolves absolute paths correctly", {
  tmp <- withr::local_tempdir()

  config_text <- paste(
    "samples:",
    "  sample1: data/s1.pod5",
    "output_dir: /absolute/path/results",
    "fasta: ref/trna.fa",
    sep = "\n"
  )
  config_path <- file.path(tmp, "config.yaml")
  writeLines(config_text, config_path)

  cfg <- read_pipeline_config(config_path)

  expect_equal(cfg$output_dir, "/absolute/path/results")
})

test_that("read_pipeline_config works with ecoli test data", {
  config_path <- clover_example("ecoli/config.yaml")
  cfg <- read_pipeline_config(config_path)

  expect_type(cfg, "list")
  expect_equal(nrow(cfg$samples), 6)
  expect_true(all(grepl("^wt-15-", cfg$samples$sample_id)))
})

test_that("list_pipeline_files constructs correct paths", {
  config <- list(
    samples = tibble::tibble(
      sample_id = c("s1", "s2"),
      data_path = c("/data/s1.pod5", "/data/s2.pod5")
    ),
    output_dir = "/out"
  )

  files <- list_pipeline_files(config, types = c("charging", "bcerror"))

  expect_s3_class(files, "tbl_df")
  expect_named(files, c("sample_id", "type", "path"))
  expect_equal(nrow(files), 4) # 2 samples x 2 types

  # Check charging path format
  charging_files <- files[files$type == "charging", ]
  expect_true(all(grepl("charging\\.cpm\\.tsv\\.gz$", charging_files$path)))

  # Check bcerror path format
  bcerror_files <- files[files$type == "bcerror", ]
  expect_true(all(grepl("bcerror\\.tsv\\.gz$", bcerror_files$path)))
})

test_that("list_pipeline_files validates types", {
  config <- list(
    samples = tibble::tibble(sample_id = "s1", data_path = "/d/s1"),
    output_dir = "/out"
  )

  expect_error(
    list_pipeline_files(config, types = "invalid_type"),
    "should be one of"
  )
})
