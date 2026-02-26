# Odds ratio downstream analysis -----------------------------------------------

#' Clean odds ratio data.
#'
#' Prepare odds ratio data for downstream analysis by capping infinite
#' log odds ratio values and optionally removing adapter positions.
#'
#' @param data A tibble of odds ratio data with columns `ref`,
#'   `odds_ratio`, and `log_odds_ratio`.
#' @param cap_inf Numeric value added to (or subtracted from) the
#'   finite range when capping positive (or negative) infinite log
#'   odds ratios. Default `2`.
#'
#' @return A tibble with a new `log_or_clean` column where infinite
#'   values have been capped.
#'
#' @export
#'
#' @examples
#' df <- tibble::tibble(
#'   ref = "tRNA-Ala-AGC-1-1",
#'   pos1 = c(20, 34),
#'   pos2 = c(34, 58),
#'   odds_ratio = c(2.5, Inf),
#'   log_odds_ratio = c(0.92, Inf),
#'   p_value = c(0.001, 0.01),
#'   total_obs = c(200, 150)
#' )
#' clean_odds_ratios(df)
clean_odds_ratios <- function(data, cap_inf = 2) {
  finite_vals <- data$log_odds_ratio[is.finite(data$log_odds_ratio)]

  if (length(finite_vals) == 0) {
    data$log_or_clean <- data$log_odds_ratio
    return(data)
  }

  upper <- max(finite_vals, na.rm = TRUE) + cap_inf
  lower <- min(finite_vals, na.rm = TRUE) - cap_inf

  dplyr::mutate(
    data,
    log_or_clean = dplyr::case_when(
      is.infinite(.data$log_odds_ratio) & .data$log_odds_ratio > 0 ~ upper,
      is.infinite(.data$log_odds_ratio) & .data$log_odds_ratio < 0 ~ lower,
      .default = .data$log_odds_ratio
    )
  )
}

#' Filter odds ratios for structure linkage arcs.
#'
#' Convenience filter for odds ratio data that returns a tibble ready
#' for the `linkages` parameter of [plot_tRNA_structure()]. Filters
#' by BH-adjusted p-value, observation count, and log odds ratio
#' magnitude, then selects the columns needed for plotting.
#'
#' @param data A tibble of odds ratio data, typically from
#'   [clean_odds_ratios()], with columns `pos1`, `pos2`,
#'   `log_odds_ratio`, `p_adjusted`, and `total_obs`.
#' @param max_p Maximum adjusted p-value (BH FDR) to retain. Default `0.01`.
#' @param min_obs Minimum total observations to retain. Default `100`.
#' @param min_lor Minimum absolute log odds ratio to retain. Default
#'   `1.0`.
#'
#' @return A tibble with columns `pos1`, `pos2`, and `value` (the log
#'   odds ratio), ready for [plot_tRNA_structure()].
#'
#' @export
#'
#' @examples
#' df <- tibble::tibble(
#'   pos1 = c(20, 34, 10),
#'   pos2 = c(34, 58, 45),
#'   odds_ratio = c(4.0, 0.3, 1.1),
#'   log_odds_ratio = c(1.4, -1.2, 0.1),
#'   p_adjusted = c(0.001, 0.005, 0.5),
#'   total_obs = c(200, 150, 50)
#' )
#' filter_linkages(df)
filter_linkages <- function(data, max_p = 0.01, min_obs = 100, min_lor = 1.0) {
  data |>
    dplyr::filter(
      .data$p_adjusted < max_p,
      .data$total_obs >= min_obs,
      abs(.data$log_odds_ratio) >= min_lor
    ) |>
    dplyr::transmute(
      pos1 = .data$pos1,
      pos2 = .data$pos2,
      value = .data$log_odds_ratio
    )
}

#' Aggregate odds ratios to isodecoder level.
#'
#' Collapse per-gene odds ratios to isodecoder level by removing gene
#' copy numbers from the reference name and summarizing across copies.
#'
#' @param data A tibble with columns `ref`, `pos1`, `pos2`,
#'   `log_or_clean` (from [clean_odds_ratios()]), `p_value`, and
#'   `total_obs`.
#' @param pattern Regular expression used to strip gene copy numbers
#'   from `ref`. Default `"-\\d+-\\d+$"` matches trailing
#'   `-<copy>-<gene>` suffixes.
#'
#' @return A tibble with columns: `isodecoder`, `pos1`, `pos2`,
#'   `mean_or`, `mean_log_or`, `sd_log_or`, `min_pval`, `total_reads`,
#'   and `n_copies`.
#'
#' @export
#'
#' @examples
#' df <- tibble::tibble(
#'   ref = c("tRNA-Ala-AGC-1-1", "tRNA-Ala-AGC-2-1"),
#'   pos1 = c(20, 20),
#'   pos2 = c(34, 34),
#'   odds_ratio = c(2.5, 3.0),
#'   log_or_clean = c(0.92, 1.10),
#'   p_value = c(0.001, 0.01),
#'   total_obs = c(200, 150)
#' )
#' aggregate_or_isodecoder(df)
aggregate_or_isodecoder <- function(data, pattern = "-\\d+-\\d+$") {
  data |>
    dplyr::mutate(
      isodecoder = stringr::str_replace(.data$ref, pattern, "")
    ) |>
    dplyr::group_by(.data$isodecoder, .data$pos1, .data$pos2) |>
    dplyr::summarise(
      mean_or = mean(.data$odds_ratio, na.rm = TRUE),
      mean_log_or = mean(.data$log_or_clean, na.rm = TRUE),
      sd_log_or = stats::sd(.data$log_or_clean, na.rm = TRUE),
      min_pval = min(.data$p_value, na.rm = TRUE),
      total_reads = sum(.data$total_obs, na.rm = TRUE),
      n_copies = dplyr::n(),
      .groups = "drop"
    )
}

#' Compute relative odds ratio at isodecoder level.
#'
#' Compare isodecoder-level odds ratios between two conditions by
#' computing the relative odds ratio (ROR) with z-score significance
#' testing.
#'
#' @param data_num A tibble of isodecoder-aggregated odds ratios for
#'   the numerator condition (from [aggregate_or_isodecoder()]).
#' @param data_den A tibble of isodecoder-aggregated odds ratios for
#'   the denominator condition.
#' @param ror_cap Maximum absolute value for the ROR. Values beyond
#'   this are capped. Default `10`.
#' @param p_method Method for p-value adjustment, passed to
#'   [stats::p.adjust()]. Default `"BH"`.
#' @param alpha Significance level for the `significant` column.
#'   Default `0.05`.
#'
#' @return A tibble with columns: `isodecoder`, `pos1`, `pos2`,
#'   `mean_log_or_num`, `mean_log_or_den`, `se_num`, `se_den`, `ror`,
#'   `ror_se`, `z_score`, `p_value`, `p_adj`, `ci_lower`, `ci_upper`,
#'   and `significant`.
#'
#' @export
#'
#' @examples
#' num <- tibble::tibble(
#'   isodecoder = rep("tRNA-Ala-AGC", 2),
#'   pos1 = c(20, 34), pos2 = c(34, 58),
#'   mean_or = c(3.0, 1.5), mean_log_or = c(1.1, 0.4),
#'   sd_log_or = c(0.2, 0.3), min_pval = c(0.001, 0.01),
#'   total_reads = c(500, 300), n_copies = c(2, 2)
#' )
#' den <- tibble::tibble(
#'   isodecoder = rep("tRNA-Ala-AGC", 2),
#'   pos1 = c(20, 34), pos2 = c(34, 58),
#'   mean_or = c(1.5, 1.2), mean_log_or = c(0.4, 0.18),
#'   sd_log_or = c(0.15, 0.25), min_pval = c(0.01, 0.05),
#'   total_reads = c(400, 250), n_copies = c(2, 2)
#' )
#' compute_ror_isodecoder(num, den)
compute_ror_isodecoder <- function(
  data_num,
  data_den,
  ror_cap = 10,
  p_method = "BH",
  alpha = 0.05
) {
  joined <- dplyr::inner_join(
    dplyr::select(
      data_num,
      "isodecoder",
      "pos1",
      "pos2",
      "mean_log_or_num" = "mean_log_or",
      "se_num" = "sd_log_or"
    ),
    dplyr::select(
      data_den,
      "isodecoder",
      "pos1",
      "pos2",
      "mean_log_or_den" = "mean_log_or",
      "se_den" = "sd_log_or"
    ),
    by = c("isodecoder", "pos1", "pos2")
  )

  joined |>
    dplyr::mutate(
      ror = .data$mean_log_or_num - .data$mean_log_or_den,
      ror = dplyr::case_when(
        .data$ror > ror_cap ~ ror_cap,
        .data$ror < -ror_cap ~ -ror_cap,
        .default = .data$ror
      ),
      ror_se = sqrt(.data$se_num^2 + .data$se_den^2),
      z_score = .data$ror / .data$ror_se,
      p_value = 2 * stats::pnorm(-abs(.data$z_score)),
      p_adj = stats::p.adjust(.data$p_value, method = p_method),
      ci_lower = .data$ror - 1.96 * .data$ror_se,
      ci_upper = .data$ror + 1.96 * .data$ror_se,
      significant = .data$p_adj < alpha
    )
}

#' Compute ratio of odds ratios between conditions.
#'
#' Compare modification co-occurrence between two conditions by computing
#' the ratio of odds ratios (ROR). Replicates within each condition are
#' aggregated using the specified function.
#'
#' @param odds_data A combined tibble of odds ratio data with a `condition`
#'   column (or column specified by `condition_col`) and `sample_id`.
#' @param condition_col Column name (string) for condition labels.
#'   Default `"condition"`.
#' @param numerator Value of `condition_col` for the numerator condition.
#' @param denominator Value of `condition_col` for the denominator condition.
#' @param min_obs Minimum `total_obs` for a pair to be included.
#'   Default `100`.
#' @param agg_fun Function to aggregate replicate log odds ratios.
#'   Default `mean`.
#'
#' @return A tibble with columns: `pos1`, `pos2`, `or_numerator`,
#'   `or_denominator`, `ror`, and `log_ror`.
#'
#' @export
#'
#' @examples
#' results <- read_pipeline_results(
#'   clover_example("ecoli/config.yaml"),
#'   types = "odds_ratios"
#' )
#' or_data <- results$odds_ratios
#' or_data$condition <- ifelse(
#'   grepl("ctl", or_data$sample_id), "ctl", "inf"
#' )
#' compute_ror(or_data, numerator = "inf", denominator = "ctl")
compute_ror <- function(
  odds_data,
  condition_col = "condition",
  numerator,
  denominator,
  min_obs = 100,
  agg_fun = mean
) {
  # Filter by minimum observations
  filtered <- odds_data |>
    dplyr::filter(.data$total_obs >= min_obs)

  # Aggregate replicates within each condition
  agg <- filtered |>
    dplyr::group_by(
      .data[[condition_col]],
      .data$pos1,
      .data$pos2
    ) |>
    dplyr::summarise(
      mean_log_or = agg_fun(.data$log_odds_ratio),
      .groups = "drop"
    )

  # Separate numerator and denominator
  num <- agg |>
    dplyr::filter(.data[[condition_col]] == numerator) |>
    dplyr::select("pos1", "pos2", "or_numerator" = "mean_log_or")

  denom <- agg |>
    dplyr::filter(.data[[condition_col]] == denominator) |>
    dplyr::select("pos1", "pos2", "or_denominator" = "mean_log_or")

  # Join and compute ROR
  dplyr::inner_join(num, denom, by = c("pos1", "pos2")) |>
    dplyr::mutate(
      log_ror = .data$or_numerator - .data$or_denominator,
      ror = exp(.data$log_ror)
    )
}

#' Compute pairwise modification co-occurrence odds ratios.
#'
#' Given a modkit `mod_calls.tsv.gz` file, build a per-read binary
#' modification matrix for each tRNA and compute Fisher's exact test
#' for each pair of positions. This is computationally expensive for
#' full datasets; use the `refs` parameter to restrict to specific
#' tRNAs.
#'
#' @param mod_calls_path Path to a `{sample}.mod_calls.tsv.gz` file
#'   from modkit.
#' @param refs Optional character vector of reference names to
#'   include. If `NULL`, all references are processed.
#' @param min_reads Minimum number of reads required for a tRNA to be
#'   included. Default `10`.
#'
#' @return A tibble with columns: `ref`, `pos1`, `pos2`,
#'   `odds_ratio`, `log_odds_ratio`, `p_value`, `total_obs`.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' compute_odds_ratios("sample.mod_calls.tsv.gz")
#' }
compute_odds_ratios <- function(mod_calls_path, refs = NULL, min_reads = 10) {
  mc <- readr::read_tsv(mod_calls_path, show_col_types = FALSE) |>
    dplyr::filter(.data$within_alignment == TRUE)

  if (!is.null(refs)) {
    mc <- mc |> dplyr::filter(.data$chrom %in% refs)
  }

  all_refs <- unique(mc$chrom)
  results <- list()

  for (r in all_refs) {
    ref_data <- mc |> dplyr::filter(.data$chrom == r)

    # Build read x position binary matrix (modified = call_code != "-")
    mat_data <- ref_data |>
      dplyr::mutate(modified = as.integer(.data$call_code != "-")) |>
      dplyr::select("read_id", "ref_position", "modified") |>
      dplyr::group_by(.data$read_id, .data$ref_position) |>
      dplyr::summarize(modified = max(.data$modified), .groups = "drop") |>
      tidyr::pivot_wider(
        names_from = "ref_position",
        values_from = "modified",
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
