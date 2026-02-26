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
      is.infinite(log_odds_ratio) & log_odds_ratio > 0 ~ upper,
      is.infinite(log_odds_ratio) & log_odds_ratio < 0 ~ lower,
      .default = log_odds_ratio
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
      isodecoder = stringr::str_replace(ref, pattern, "")
    ) |>
    dplyr::group_by(isodecoder, pos1, pos2) |>
    dplyr::summarise(
      mean_or = mean(odds_ratio, na.rm = TRUE),
      mean_log_or = mean(log_or_clean, na.rm = TRUE),
      sd_log_or = stats::sd(log_or_clean, na.rm = TRUE),
      min_pval = min(p_value, na.rm = TRUE),
      total_reads = sum(total_obs, na.rm = TRUE),
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
      isodecoder,
      pos1,
      pos2,
      mean_log_or_num = mean_log_or,
      se_num = sd_log_or
    ),
    dplyr::select(
      data_den,
      isodecoder,
      pos1,
      pos2,
      mean_log_or_den = mean_log_or,
      se_den = sd_log_or
    ),
    by = c("isodecoder", "pos1", "pos2")
  )

  joined |>
    dplyr::mutate(
      ror = mean_log_or_num - mean_log_or_den,
      ror = dplyr::case_when(
        ror > ror_cap ~ ror_cap,
        ror < -ror_cap ~ -ror_cap,
        .default = ror
      ),
      ror_se = sqrt(se_num^2 + se_den^2),
      z_score = ror / ror_se,
      p_value = 2 * stats::pnorm(-abs(z_score)),
      p_adj = stats::p.adjust(p_value, method = p_method),
      ci_lower = ror - 1.96 * ror_se,
      ci_upper = ror + 1.96 * ror_se,
      significant = p_adj < alpha
    )
}
