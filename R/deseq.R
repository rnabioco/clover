# DESeq2 wrapper functions ----------------------------------------------------

#' Build a tRNA abundance count matrix from charging data.
#'
#' Pivot charging data into an integer matrix suitable for DESeq2, where
#' each row is a tRNA and each column is a sample. Values are total
#' abundance (charged + uncharged counts).
#'
#' @param charging_data A tibble of combined charging data with columns
#'   `tRNA`, `counts_charged`, `counts_uncharged`, and `sample_id`
#'   (as returned by [read_charging_multi()]).
#' @param min_count Minimum total count across all samples for a tRNA
#'   to be retained. Default `10`.
#'
#' @return An integer matrix with tRNA names as row names and sample IDs
#'   as column names.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' charging <- read_charging_multi(paths)
#' mat <- abundance_count_matrix(charging)
#' }
abundance_count_matrix <- function(charging_data, min_count = 10) {
  abundance <- charging_data |>
    dplyr::mutate(
      total_count = as.integer(
        round(counts_charged + counts_uncharged)
      )
    ) |>
    dplyr::select(tRNA, sample_id, total_count) |>
    tidyr::pivot_wider(
      names_from = sample_id,
      values_from = total_count,
      values_fill = 0L
    )

  mat <- as.matrix(abundance[, -1, drop = FALSE])
  rownames(mat) <- abundance$tRNA
  storage.mode(mat) <- "integer"

  # Filter low-count tRNAs
  keep <- rowSums(mat) >= min_count
  mat[keep, , drop = FALSE]
}

#' Build a charging count matrix for differential charging analysis.
#'
#' Create a count matrix with two columns per sample (charged and uncharged),
#' suitable for a DESeq2 interaction design that tests for differences in
#' charging ratios between conditions.
#'
#' @param charging_data A tibble of combined charging data with columns
#'   `tRNA`, `counts_charged`, `counts_uncharged`, and `sample_id`
#'   (as returned by [read_charging_multi()]).
#' @param min_count Minimum total count across all columns for a tRNA
#'   to be retained. Default `10`.
#'
#' @return An integer matrix with tRNA names as row names. Column names
#'   are `{sample_id}_charged` and `{sample_id}_uncharged`.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' charging <- read_charging_multi(paths)
#' mat <- charging_count_matrix(charging)
#' }
charging_count_matrix <- function(charging_data, min_count = 10) {
  long <- charging_data |>
    dplyr::mutate(
      counts_charged = as.integer(round(counts_charged)),
      counts_uncharged = as.integer(round(counts_uncharged))
    ) |>
    tidyr::pivot_longer(
      cols = c(counts_charged, counts_uncharged),
      names_to = "charge_status",
      values_to = "count"
    ) |>
    dplyr::mutate(
      charge_status = sub("^counts_", "", charge_status),
      col_name = paste0(sample_id, "_", charge_status)
    ) |>
    dplyr::select(tRNA, col_name, count) |>
    tidyr::pivot_wider(
      names_from = col_name,
      values_from = count,
      values_fill = 0L
    )

  mat <- as.matrix(long[, -1, drop = FALSE])
  rownames(mat) <- long$tRNA
  storage.mode(mat) <- "integer"

  # Filter low-count tRNAs
  keep <- rowSums(mat) >= min_count
  mat[keep, , drop = FALSE]
}

#' Build column metadata for DESeq2.
#'
#' Create a `data.frame` of sample metadata suitable for
#' [DESeq2::DESeqDataSetFromMatrix()]. For charging count matrices,
#' the `charge_status` factor is added automatically.
#'
#' @param count_matrix A count matrix from [abundance_count_matrix()] or
#'   [charging_count_matrix()].
#' @param sample_info An optional data frame with a `sample_id` column and
#'   additional experimental factor columns (e.g., `condition`, `replicate`).
#'   If `NULL`, a minimal data frame is created from column names.
#'
#' @return A data frame with row names matching `colnames(count_matrix)`.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' mat <- abundance_count_matrix(charging)
#' sample_info <- data.frame(
#'   sample_id = c("wt_1", "wt_2", "mut_1", "mut_2"),
#'   condition = c("wt", "wt", "mut", "mut")
#' )
#' coldata <- build_coldata(mat, sample_info)
#' }
build_coldata <- function(count_matrix, sample_info = NULL) {
  col_names <- colnames(count_matrix)

  # Detect charging matrix format
  is_charging <- any(grepl("_(charged|uncharged)$", col_names))

  if (is_charging) {
    # Extract sample_id and charge_status from column names
    coldata <- data.frame(
      sample_id = sub("_(charged|uncharged)$", "", col_names),
      charge_status = factor(
        sub("^.*_(charged|uncharged)$", "\\1", col_names),
        levels = c("uncharged", "charged")
      ),
      row.names = col_names,
      stringsAsFactors = FALSE
    )
  } else {
    coldata <- data.frame(
      sample_id = col_names,
      row.names = col_names,
      stringsAsFactors = FALSE
    )
  }

  # Merge user-supplied sample info
  if (!is.null(sample_info)) {
    sample_info <- as.data.frame(sample_info)
    coldata <- merge(
      coldata, sample_info,
      by = "sample_id", all.x = TRUE, sort = FALSE
    )
    rownames(coldata) <- col_names
  }

  coldata
}

#' Run DESeq2 differential analysis.
#'
#' Wrapper around [DESeq2::DESeqDataSetFromMatrix()] and [DESeq2::DESeq()]
#' for tRNA abundance or charging data.
#'
#' @param count_matrix An integer count matrix (from [abundance_count_matrix()]
#'   or [charging_count_matrix()]).
#' @param coldata A data frame of sample metadata (from [build_coldata()]).
#' @param design A formula specifying the design (e.g., `~ condition` or
#'   `~ condition + charge_status + condition:charge_status`).
#' @param ... Additional arguments passed to [DESeq2::DESeq()].
#'
#' @return A `DESeqDataSet` object.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' dds <- run_deseq(mat, coldata, design = ~ condition)
#' }
run_deseq <- function(count_matrix, coldata, design, ...) {
  rlang::check_installed("DESeq2", reason = "to run differential analysis.")

  dds <- DESeq2::DESeqDataSetFromMatrix(
    countData = count_matrix,
    colData = coldata,
    design = design
  )

  DESeq2::DESeq(dds, ...)
}

#' Tidy DESeq2 results into a tibble.
#'
#' Extract results from a DESeq2 analysis and return a tidy tibble with
#' tRNA identifiers and significance flags.
#'
#' @param dds A `DESeqDataSet` object (from [run_deseq()]).
#' @param contrast A contrast specification: either a character vector of
#'   length 3 (e.g., `c("condition", "mutant", "wildtype")`) or a list
#'   for coefficient-based contrasts.
#' @param padj_cutoff Adjusted p-value threshold for significance.
#'   Default `0.05`.
#'
#' @return A tibble with columns: `tRNA`, `log2FoldChange`, `lfcSE`,
#'   `pvalue`, `padj`, and `significant` (logical).
#'
#' @export
#'
#' @examples
#' \dontrun{
#' res <- tidy_deseq_results(dds, contrast = c("condition", "mut", "wt"))
#' res
#' }
tidy_deseq_results <- function(dds, contrast, padj_cutoff = 0.05) {
  rlang::check_installed("DESeq2", reason = "to extract results.")

  res <- DESeq2::results(dds, contrast = contrast)
  res_df <- as.data.frame(res)

  tibble::tibble(
    tRNA = rownames(res_df),
    log2FoldChange = res_df$log2FoldChange,
    lfcSE = res_df$lfcSE,
    pvalue = res_df$pvalue,
    padj = res_df$padj,
    significant = !is.na(res_df$padj) & res_df$padj < padj_cutoff
  )
}
