# Charging and odds ratio reading functions ------------------------------------

#' Read a charging CPM file.
#'
#' Read a per-tRNA charging CPM file produced by the tRNA sequencing pipeline.
#' These files contain counts and CPM values for charged and uncharged reads.
#'
#' @param path Path to a `{sample}.charging.cpm.tsv.gz` file.
#'
#' @return A tibble with columns including `ref`, `counts_charged`,
#'   `counts_uncharged`, `cpm_charged`, `cpm_uncharged`, and `total_count`.
#'
#' @export
#'
#' @examples
#' path <- clover_example(
#'   "ecoli/summary/tables/wt-15-ctl-01/wt-15-ctl-01.charging.cpm.tsv.gz"
#' )
#' read_charging(path)
read_charging <- function(path) {
  readr::read_tsv(path, show_col_types = FALSE) |>
    dplyr::rename("ref" = "tRNA")
}

#' Read an odds ratios file.
#'
#' Read a per-position-pair odds ratio file produced by the tRNA sequencing
#' pipeline. These files contain pairwise modification co-occurrence statistics.
#'
#' @param path Path to a `{sample}.odds_ratios_filtered.tsv.gz` file.
#'
#' @return A tibble with columns including `ref`, `pos1`, `pos2`,
#'   `odds_ratio`, `log_odds_ratio`, `p_value`, and `total_obs`.
#'
#' @export
#'
#' @examples
#' path <- clover_example(
#'   "ecoli/summary/tables/wt-15-ctl-01/wt-15-ctl-01.odds_ratios_filtered.tsv.gz"
#' )
#' read_odds_ratios(path)
read_odds_ratios <- function(path) {
  col_types <- readr::cols(
    tRNA = readr::col_character(),
    pos1 = readr::col_character(),
    pos2 = readr::col_character(),
    n00 = readr::col_double(),
    n01 = readr::col_double(),
    n10 = readr::col_double(),
    n11 = readr::col_double(),
    total_obs = readr::col_double(),
    odds_ratio = readr::col_double(),
    log_odds_ratio = readr::col_double(),
    se_log_or = readr::col_double(),
    ci_lower = readr::col_double(),
    ci_upper = readr::col_double(),
    fisher_or = readr::col_double(),
    p_value = readr::col_double(),
    p_adjusted = readr::col_double()
  )
  d <- readr::read_tsv(path, col_types = col_types)
  if ("tRNA" %in% names(d)) {
    d <- dplyr::rename(d, "ref" = "tRNA")
  }
  d
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
#' paths <- c(
#'   "wt-15-ctl-01" = clover_example(
#'     "ecoli/summary/tables/wt-15-ctl-01/wt-15-ctl-01.charging.cpm.tsv.gz"
#'   ),
#'   "wt-15-ctl-02" = clover_example(
#'     "ecoli/summary/tables/wt-15-ctl-02/wt-15-ctl-02.charging.cpm.tsv.gz"
#'   )
#' )
#' read_charging_multi(paths)
read_charging_multi <- function(paths) {
  read_multi(paths, read_charging)
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
#' or_file1 <- paste0(
#'   "ecoli/summary/tables/wt-15-ctl-01/",
#'   "wt-15-ctl-01.odds_ratios_filtered.tsv.gz"
#' )
#' or_file2 <- paste0(
#'   "ecoli/summary/tables/wt-15-ctl-02/",
#'   "wt-15-ctl-02.odds_ratios_filtered.tsv.gz"
#' )
#' paths <- c(
#'   "wt-15-ctl-01" = clover_example(or_file1),
#'   "wt-15-ctl-02" = clover_example(or_file2)
#' )
#' read_odds_ratios_multi(paths)
read_odds_ratios_multi <- function(paths) {
  read_multi(paths, read_odds_ratios)
}

# Internal helpers -----------------------------------------------------------

read_multi <- function(paths, reader) {
  if (is.null(names(paths))) {
    cli_abort("{.arg paths} must be a named character vector.")
  }

  purrr::map(paths, reader) |>
    dplyr::bind_rows(.id = "sample_id")
}

# Charging ratio comparison ------------------------------------------------------

#' Compute charging ratio differences between conditions.
#'
#' Compare per-tRNA charging ratios (charged / total) between two
#' conditions, summarizing replicates with mean and standard error and
#' computing the between-condition difference with propagated SE.
#'
#' @param charging_data A tibble from [read_charging_multi()] with an
#'   added condition column. Must contain `ref`, `counts_charged`,
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
#'   \item{ref}{tRNA identifier (factor ordered by `diff`).}
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
#' config_path <- clover_example("ecoli/config.yaml")
#' results <- read_pipeline_results(config_path, types = "charging")
#' charging <- results$charging
#' charging$condition <- ifelse(
#'   grepl("ctl", charging$sample_id), "ctl", "inf"
#' )
#' compute_charging_diffs(
#'   charging,
#'   numerator = "inf",
#'   denominator = "ctl"
#' )
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
    dplyr::filter(!grepl("-uncharged$", .data$ref)) |>
    dplyr::mutate(
      total = .data$counts_charged + .data$counts_uncharged,
      charging_ratio = .data$counts_charged / .data$total
    ) |>
    dplyr::filter(.data$total >= min_count)

  # Optionally keep only top N tRNAs by total abundance
  if (!is.null(n_top)) {
    top_refs <- ratios |>
      dplyr::group_by(.data$ref) |>
      dplyr::summarise(total = sum(.data$total), .groups = "drop") |>
      dplyr::slice_max(.data$total, n = n_top) |>
      dplyr::pull(.data$ref)

    ratios <- ratios |>
      dplyr::filter(.data$ref %in% top_refs)
  }

  # Summarize by tRNA and condition
  ratio_summary <- ratios |>
    dplyr::group_by(.data$ref, .data[[condition_col]]) |>
    dplyr::summarise(
      mean_ratio = mean(.data$charging_ratio),
      se_ratio = stats::sd(.data$charging_ratio) / sqrt(dplyr::n()),
      .groups = "drop"
    )

  # Pivot wider by condition
  wide <- ratio_summary |>
    tidyr::pivot_wider(
      names_from = dplyr::all_of(condition_col),
      values_from = c("mean_ratio", "se_ratio")
    )

  num_ratio <- paste0("mean_ratio_", numerator)
  den_ratio <- paste0("mean_ratio_", denominator)
  num_se <- paste0("se_ratio_", numerator)
  den_se <- paste0("se_ratio_", denominator)

  wide |>
    dplyr::transmute(
      ref = .data$ref,
      ratio_numerator = .data[[num_ratio]],
      ratio_denominator = .data[[den_ratio]],
      se_numerator = .data[[num_se]],
      se_denominator = .data[[den_se]],
      diff = .data[[num_ratio]] - .data[[den_ratio]],
      se_diff = sqrt(.data[[num_se]]^2 + .data[[den_se]]^2)
    ) |>
    dplyr::filter(!is.na(.data$diff)) |>
    dplyr::mutate(ref = forcats::fct_reorder(.data$ref, .data$diff))
}
