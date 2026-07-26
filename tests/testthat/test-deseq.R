test_that("abundance_count_matrix creates correct matrix", {
  charging <- tibble::tibble(
    ref = rep(c("tRNA-Ala", "tRNA-Gly", "tRNA-Ser"), 2),
    counts_charged = c(100, 200, 5, 120, 180, 3),
    counts_uncharged = c(50, 80, 2, 60, 90, 1),
    sample_id = rep(c("s1", "s2"), each = 3)
  )

  mat <- abundance_count_matrix(charging, min_count = 10)

  expect_true(is.matrix(mat))
  expect_true(is.integer(mat))
  expect_equal(ncol(mat), 2)
  expect_equal(colnames(mat), c("s1", "s2"))
  # tRNA-Ser has total 5+2+3+1=11 >= 10, so it should be included
  expect_equal(nrow(mat), 3)
  # Check values: total = charged + uncharged
  expect_equal(mat["tRNA-Ala", "s1"], 150L)
  expect_equal(mat["tRNA-Gly", "s2"], 270L)
})

test_that("abundance_count_matrix filters low-count tRNAs", {
  charging <- tibble::tibble(
    ref = rep(c("tRNA-Ala", "tRNA-Low"), 2),
    counts_charged = c(100, 1, 120, 2),
    counts_uncharged = c(50, 0, 60, 1),
    sample_id = rep(c("s1", "s2"), each = 2)
  )

  mat <- abundance_count_matrix(charging, min_count = 10)

  expect_equal(nrow(mat), 1)
  expect_equal(rownames(mat), "tRNA-Ala")
})

test_that("charging_count_matrix creates correct matrix", {
  charging <- tibble::tibble(
    ref = rep(c("tRNA-Ala", "tRNA-Gly"), 2),
    counts_charged = c(100, 200, 120, 180),
    counts_uncharged = c(50, 80, 60, 90),
    sample_id = rep(c("s1", "s2"), each = 2)
  )

  mat <- charging_count_matrix(charging, min_count = 10)

  expect_true(is.matrix(mat))
  expect_true(is.integer(mat))
  expect_equal(nrow(mat), 2)
  # Should have 4 columns: s1_charged, s1_uncharged, s2_charged, s2_uncharged
  expect_equal(ncol(mat), 4)
  expect_true(all(c("s1_charged", "s1_uncharged") %in% colnames(mat)))
})

test_that("build_coldata works for abundance matrix", {
  mat <- matrix(1:6, nrow = 2, ncol = 3)
  colnames(mat) <- c("s1", "s2", "s3")

  coldata <- build_coldata(mat)

  expect_s3_class(coldata, "data.frame")
  expect_equal(rownames(coldata), c("s1", "s2", "s3"))
  expect_true("sample_id" %in% names(coldata))
})

test_that("build_coldata detects charging matrix format", {
  mat <- matrix(1:8, nrow = 2, ncol = 4)
  colnames(mat) <- c("s1_charged", "s1_uncharged", "s2_charged", "s2_uncharged")

  coldata <- build_coldata(mat)

  expect_true("charge_status" %in% names(coldata))
  expect_s3_class(coldata$charge_status, "factor")
  expect_equal(levels(coldata$charge_status), c("uncharged", "charged"))
})

test_that("build_coldata merges sample info", {
  mat <- matrix(1:4, nrow = 2, ncol = 2)
  colnames(mat) <- c("s1", "s2")

  sample_info <- data.frame(
    sample_id = c("s1", "s2"),
    condition = c("wt", "mut")
  )

  coldata <- build_coldata(mat, sample_info)

  expect_true("condition" %in% names(coldata))
  expect_equal(coldata$condition, c("wt", "mut"))
})

test_that("run_deseq and tidy_deseq_results work end-to-end", {
  skip_if_not_installed("DESeq2")

  # Create synthetic count data with enough genes for dispersion estimation
  set.seed(42)
  n_genes <- 20
  n_samples <- 6
  lambdas <- sample(50:500, n_genes, replace = TRUE)

  counts <- matrix(
    as.integer(rpois(n_genes * n_samples, lambda = rep(lambdas, n_samples))),
    nrow = n_genes,
    ncol = n_samples
  )
  rownames(counts) <- paste0("tRNA-", seq_len(n_genes))
  colnames(counts) <- c("wt1", "wt2", "wt3", "mut1", "mut2", "mut3")

  coldata <- data.frame(
    sample_id = colnames(counts),
    condition = factor(rep(c("wt", "mut"), each = 3)),
    row.names = colnames(counts)
  )

  # Suppress DESeq2 dispersion warnings on small synthetic datasets
  suppressWarnings(
    dds <- run_deseq(counts, coldata, design = ~condition)
  )
  expect_s4_class(dds, "DESeqDataSet")

  res <- tidy_deseq_results(dds, contrast = c("condition", "mut", "wt"))
  expect_s3_class(res, "tbl_df")
  expect_named(
    res,
    c("ref", "log2FoldChange", "lfcSE", "pvalue", "padj", "significant")
  )
  expect_type(res$significant, "logical")
  expect_equal(nrow(res), n_genes)
})

# A small fitted model to extract coefficients from. With `interaction = TRUE`
# the columns are split into charged and uncharged halves, which is the shape of
# a differential charging design.
synthetic_dds <- function(interaction = FALSE, n_genes = 20) {
  set.seed(42)
  n_samples <- 6
  lambdas <- sample(50:500, n_genes, replace = TRUE)

  counts <- matrix(
    as.integer(rpois(n_genes * n_samples, lambda = rep(lambdas, n_samples))),
    nrow = n_genes,
    ncol = n_samples
  )
  rownames(counts) <- paste0("tRNA-", seq_len(n_genes))
  colnames(counts) <- c("wt1", "wt2", "wt3", "mut1", "mut2", "mut3")

  coldata <- data.frame(
    sample_id = colnames(counts),
    condition = factor(rep(c("wt", "mut"), each = 3), levels = c("wt", "mut")),
    row.names = colnames(counts)
  )

  if (!interaction) {
    return(suppressWarnings(run_deseq(counts, coldata, design = ~condition)))
  }

  counts <- cbind(counts, counts)
  colnames(counts) <- c(
    paste0(colnames(coldata)[0], rownames(coldata), "_uncharged"),
    paste0(rownames(coldata), "_charged")
  )

  coldata <- data.frame(
    condition = factor(rep(coldata$condition, 2), levels = c("wt", "mut")),
    charge_status = factor(
      rep(c("uncharged", "charged"), each = n_samples),
      levels = c("uncharged", "charged")
    ),
    row.names = colnames(counts)
  )

  suppressWarnings(
    run_deseq(counts, coldata, design = ~ condition * charge_status)
  )
}

test_that("tidy_deseq_results extracts a coefficient by name", {
  skip_if_not_installed("DESeq2")

  dds <- synthetic_dds()
  coef <- setdiff(DESeq2::resultsNames(dds), "Intercept")[1]

  by_name <- tidy_deseq_results(dds, name = coef)

  expect_s3_class(by_name, "tbl_df")
  expect_named(
    by_name,
    c("ref", "log2FoldChange", "lfcSE", "pvalue", "padj", "significant")
  )
  # The single coefficient of a two-level factor is that factor's contrast, so
  # the two routes have to agree.
  by_contrast <- tidy_deseq_results(dds, contrast = c("condition", "mut", "wt"))
  expect_equal(by_name$log2FoldChange, by_contrast$log2FoldChange)
  expect_equal(by_name$pvalue, by_contrast$pvalue)
})

test_that("tidy_deseq_results reaches an interaction coefficient", {
  skip_if_not_installed("DESeq2")

  # An interaction term is unreachable through `contrast`, which is the reason
  # `name` exists: this is the shape of a differential charging model.
  dds <- synthetic_dds(interaction = TRUE)
  coef <- grep("\\.", DESeq2::resultsNames(dds), value = TRUE)
  expect_length(coef, 1)

  res <- tidy_deseq_results(dds, name = coef)

  expect_s3_class(res, "tbl_df")
  expect_equal(nrow(res), 20)
})

test_that("tidy_deseq_results requires exactly one of contrast and name", {
  skip_if_not_installed("DESeq2")

  dds <- synthetic_dds()

  expect_error(tidy_deseq_results(dds), "exactly one")
  expect_error(
    tidy_deseq_results(
      dds,
      contrast = c("condition", "mut", "wt"),
      name = "condition_mut_vs_wt"
    ),
    "exactly one"
  )
})

test_that("tidy_deseq_results rejects an unknown coefficient name", {
  skip_if_not_installed("DESeq2")

  dds <- synthetic_dds()

  expect_error(
    tidy_deseq_results(dds, name = "no_such_coefficient"),
    "not a coefficient"
  )
})
