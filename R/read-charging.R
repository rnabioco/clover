# Charging and odds ratio reading functions ------------------------------------

#' Read a charging CPM file.
#'
#' Read a per-tRNA charging CPM file produced by the tRNA sequencing pipeline.
#' These files contain counts and CPM values for charged and uncharged reads.
#'
#' @param path Path to a `{sample}.charging.cpm.tsv.gz` file.
#'
#' @return A tibble with columns including `tRNA`, `counts_charged`,
#'   `counts_uncharged`, `cpm_charged`, `cpm_uncharged`, and `total_count`.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' charging <- read_charging("sample1.charging.cpm.tsv.gz")
#' charging
#' }
read_charging <- function(path) {
  readr::read_tsv(path, show_col_types = FALSE)
}

#' Read an odds ratios file.
#'
#' Read a per-position-pair odds ratio file produced by the tRNA sequencing
#' pipeline. These files contain pairwise modification co-occurrence statistics.
#'
#' @param path Path to a `{sample}.odds_ratios.tsv.gz` file.
#'
#' @return A tibble with columns including `pos1`, `pos2`, `odds_ratio`,
#'   `log_odds_ratio`, `p_value`, and `total_obs`.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' or_data <- read_odds_ratios("sample1.odds_ratios.tsv.gz")
#' or_data
#' }
read_odds_ratios <- function(path) {
  readr::read_tsv(path, show_col_types = FALSE)
}

#' Read charging CPM files for multiple samples.
#'
#' Read and combine charging CPM files from multiple samples into a single
#' tibble with a `sample_id` column.
#'
#' @param paths A named character vector of file paths. Names are used as
#'   sample identifiers.
#'
#' @return A tibble with all samples combined and a `sample_id` column.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' paths <- c(sample1 = "s1.charging.cpm.tsv.gz",
#'            sample2 = "s2.charging.cpm.tsv.gz")
#' charging <- read_charging_multi(paths)
#' }
read_charging_multi <- function(paths) {
  .read_multi(paths, read_charging)
}

#' Read odds ratio files for multiple samples.
#'
#' Read and combine odds ratio files from multiple samples into a single
#' tibble with a `sample_id` column.
#'
#' @param paths A named character vector of file paths. Names are used as
#'   sample identifiers.
#'
#' @return A tibble with all samples combined and a `sample_id` column.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' paths <- c(sample1 = "s1.odds_ratios.tsv.gz",
#'            sample2 = "s2.odds_ratios.tsv.gz")
#' or_data <- read_odds_ratios_multi(paths)
#' }
read_odds_ratios_multi <- function(paths) {
  .read_multi(paths, read_odds_ratios)
}

# Internal helpers -----------------------------------------------------------

.read_multi <- function(paths, reader) {
  if (is.null(names(paths))) {
    cli_abort("{.arg paths} must be a named character vector.")
  }

  tbls <- lapply(names(paths), function(sid) {
    tbl <- reader(paths[[sid]])
    tbl$sample_id <- sid
    tbl
  })

  dplyr::bind_rows(tbls)
}

# Charging ratio comparison ------------------------------------------------------

#' Compute charging ratio differences between conditions.
#'
#' Compare per-tRNA charging ratios (charged / total) between two
#' conditions, summarizing replicates with mean and standard error and
#' computing the between-condition difference with propagated SE.
#'
#' @param charging_data A tibble from [read_charging_multi()] with an
#'   added condition column. Must contain `tRNA`, `counts_charged`,
#'   `counts_uncharged`, `sample_id`, and the column named by
#'   `condition_col`.
#' @param condition_col Column name (string) for condition labels.
#'   Default `"condition"`.
#' @param numerator Value of `condition_col` for the numerator
#'   condition.
#' @param denominator Value of `condition_col` for the denominator
#'   condition.
#' @param min_count Minimum total reads (charged + uncharged) per tRNA
#'   per sample to include. Default `50`.
#' @param n_top If non-`NULL`, keep only the top `n_top` tRNAs by
#'   total abundance across all samples. Default `NULL` (keep all).
#'
#' @return A tibble with columns:
#' \describe{
#'   \item{tRNA}{tRNA identifier (factor ordered by `diff`).}
#'   \item{ratio_numerator}{Mean charging ratio for the numerator
#'     condition.}
#'   \item{ratio_denominator}{Mean charging ratio for the denominator
#'     condition.}
#'   \item{se_numerator}{Standard error of the numerator ratio.}
#'   \item{se_denominator}{Standard error of the denominator ratio.}
#'   \item{diff}{Difference: `ratio_numerator - ratio_denominator`.}
#'   \item{se_diff}{Propagated SE:
#'     `sqrt(se_numerator^2 + se_denominator^2)`.}
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' paths <- c(
#'   ctl1 = "ctl1.charging.cpm.tsv.gz",
#'   ctl2 = "ctl2.charging.cpm.tsv.gz",
#'   inf1 = "inf1.charging.cpm.tsv.gz",
#'   inf2 = "inf2.charging.cpm.tsv.gz"
#' )
#' charging <- read_charging_multi(paths)
#' charging$condition <- ifelse(
#'   grepl("ctl", charging$sample_id), "ctl", "inf"
#' )
#' diffs <- compute_charging_diffs(
#'   charging,
#'   numerator = "inf",
#'   denominator = "ctl"
#' )
#' }
compute_charging_diffs <- function(
  charging_data,
  condition_col = "condition",
  numerator,

  denominator,
  min_count = 50,
  n_top = NULL
) {
  if (!condition_col %in% names(charging_data)) {
    cli_abort(
      "Column {.val {condition_col}} not found in {.arg charging_data}."
    )
  }

  # Filter uncharged variants and compute per-sample charging ratio
  ratios <- charging_data |>
    dplyr::filter(!grepl("-uncharged$", tRNA)) |>
    dplyr::mutate(
      total = counts_charged + counts_uncharged,
      charging_ratio = counts_charged / total
    ) |>
    dplyr::filter(total >= min_count)

  # Optionally keep only top N tRNAs by total abundance

  if (!is.null(n_top)) {
    top_trnas <- ratios |>
      dplyr::group_by(tRNA) |>
      dplyr::summarise(total = sum(total), .groups = "drop") |>
      dplyr::slice_max(total, n = n_top) |>
      dplyr::pull(tRNA)

    ratios <- ratios |>
      dplyr::filter(tRNA %in% top_trnas)
  }

  # Summarize by tRNA and condition
  ratio_summary <- ratios |>
    dplyr::group_by(tRNA, .data[[condition_col]]) |>
    dplyr::summarise(
      mean_ratio = mean(charging_ratio),
      se_ratio = stats::sd(charging_ratio) / sqrt(dplyr::n()),
      .groups = "drop"
    )

  # Pivot wider by condition
  wide <- ratio_summary |>
    tidyr::pivot_wider(
      names_from = dplyr::all_of(condition_col),
      values_from = c(mean_ratio, se_ratio)
    )

  num_ratio <- paste0("mean_ratio_", numerator)
  den_ratio <- paste0("mean_ratio_", denominator)
  num_se <- paste0("se_ratio_", numerator)
  den_se <- paste0("se_ratio_", denominator)

  wide |>
    dplyr::transmute(
      tRNA,
      ratio_numerator = .data[[num_ratio]],
      ratio_denominator = .data[[den_ratio]],
      se_numerator = .data[[num_se]],
      se_denominator = .data[[den_se]],
      diff = .data[[num_ratio]] - .data[[den_ratio]],
      se_diff = sqrt(.data[[num_se]]^2 + .data[[den_se]]^2)
    ) |>
    dplyr::filter(!is.na(diff)) |>
    dplyr::mutate(tRNA = forcats::fct_reorder(tRNA, diff))
}

# Odds ratio computation from mod_calls -----------------------------------------

#' Compute pairwise modification co-occurrence odds ratios.
#'
#' Given a modkit `mod_calls.tsv.gz` file, build a per-read binary modification
#' matrix for each tRNA and compute Fisher's exact test for each pair of
#' positions. This is computationally expensive for full datasets; use the
#' `refs` parameter to restrict to specific tRNAs.
#'
#' @param mod_calls_path Path to a `{sample}.mod_calls.tsv.gz` file from modkit.
#' @param refs Optional character vector of reference names to include.
#'   If `NULL`, all references are processed.
#' @param min_reads Minimum number of reads required for a tRNA to be
#'   included. Default `10`.
#'
#' @return A tibble with columns: `ref`, `pos1`, `pos2`, `odds_ratio`,
#'   `log_odds_ratio`, `p_value`, `total_obs`.
#'
#' @export
compute_odds_ratios <- function(mod_calls_path, refs = NULL, min_reads = 10) {
  mc <- readr::read_tsv(mod_calls_path, show_col_types = FALSE) |>
    dplyr::filter(within_alignment == TRUE)

  if (!is.null(refs)) {
    mc <- mc |> dplyr::filter(chrom %in% refs)
  }

  all_refs <- unique(mc$chrom)
  results <- list()

  for (r in all_refs) {
    ref_data <- mc |> dplyr::filter(chrom == r)

    # Build read x position binary matrix (modified = call_code != "-")
    mat_data <- ref_data |>
      dplyr::mutate(modified = as.integer(call_code != "-")) |>
      dplyr::select(read_id, ref_position, modified) |>
      dplyr::group_by(read_id, ref_position) |>
      dplyr::summarize(modified = max(modified), .groups = "drop") |>
      tidyr::pivot_wider(
        names_from = ref_position,
        values_from = modified,
        values_fill = 0L
      )

    if (nrow(mat_data) < min_reads) {
      next
    }

    pos_cols <- setdiff(names(mat_data), "read_id")
    if (length(pos_cols) < 2) {
      next
    }

    # Compute Fisher's test for each pair via C++
    int_mat <- as.matrix(mat_data[, pos_cols, drop = FALSE])
    storage.mode(int_mat) <- "integer"

    pair_results <- pairwise_fisher_exact(int_mat)

    if (nrow(pair_results) > 0) {
      results[[length(results) + 1]] <- tibble::tibble(
        ref = r,
        pos1 = pos_cols[pair_results$pos1],
        pos2 = pos_cols[pair_results$pos2],
        odds_ratio = pair_results$odds_ratio,
        log_odds_ratio = pair_results$log_odds_ratio,
        p_value = pair_results$p_value,
        total_obs = pair_results$total_obs
      )
    }
  }

  dplyr::bind_rows(results)
}
