# Charging odds ratio analysis -------------------------------------------------

#' Read per-read charging calls.
#'
#' Read a per-read charging likelihood file produced by the tRNA sequencing
#' pipeline and binarize it at `ml_threshold`. Unlike [read_charging()], which
#' returns per-tRNA aggregate counts, this returns one row per read and is the
#' input required by [compute_charging_odds_ratios()].
#'
#' @param path Path to a `{sample}.charging_prob.tsv.gz` file, with columns
#'   `read_id`, `tRNA`, and `charging_likelihood`.
#' @param ml_threshold Minimum `charging_likelihood` (the `CL` tag) for a read
#'   to be called charged. Default `200`, matching the pipeline default.
#'
#' @return A tibble with columns `read_id`, `ref`, and `charged` (0/1).
#'
#' @seealso [compute_charging_odds_ratios()], [read_charging()]
#'
#' @export
#'
#' @examples
#' path <- clover_example("ecoli/charging_calls.tsv.gz")
#' read_charging_calls(path)
read_charging_calls <- function(path, ml_threshold = 200) {
  raw_tbl <- readr::read_tsv(path, show_col_types = FALSE)

  required <- c("read_id", "charging_likelihood")
  missing <- setdiff(required, names(raw_tbl))
  if (length(missing) > 0) {
    cli_abort(
      "Charging file is missing column{?s} {.field {missing}}."
    )
  }

  ref_col <- if ("tRNA" %in% names(raw_tbl)) "tRNA" else "ref"

  raw_tbl |>
    dplyr::transmute(
      read_id = .data$read_id,
      ref = .data[[ref_col]],
      charged = as.integer(.data$charging_likelihood >= ml_threshold)
    )
}

#' Call modified sites from base-calling error.
#'
#' Select candidate modification sites from a bcerror tibble by thresholding
#' the base-calling error rate. The result is a site list suitable for the
#' `sites` argument of [compute_charging_odds_ratios()] and
#' [compute_odds_ratios()].
#'
#' This is an alternative to selecting sites from direct modification calls:
#' base-calling error is an orthogonal, model-free signal, so it can be used
#' when the modification-caller channels are not trusted for a given dataset.
#' Note that error-derived sites carry no modification identity — a site is
#' called because reads misbehave there, not because a particular modification
#' was identified.
#'
#' @param bcerror A bcerror tibble from [read_bcerror()], with columns `ref`,
#'   `pos`, `cov`, and `error_rate`.
#' @param min_error Minimum `error_rate` for a position to be called. Default
#'   `0.1`.
#' @param min_cov Minimum coverage (`cov`) for a position to be considered.
#'   Default `20`.
#' @param refs Optional character vector of reference names to restrict to.
#'   If `NULL` (default), all references are considered.
#'
#' @return A tibble with columns `ref`, `pos`, `cov`, and `error_rate`, sorted
#'   by reference and position.
#'
#' @seealso [read_bcerror()], [compute_charging_odds_ratios()]
#'
#' @export
#'
#' @examples
#' bcerr_path <- clover_example(
#'   "ecoli/summary/tables/wt-15-ctl-01/wt-15-ctl-01.bcerror.tsv.gz"
#' )
#' bcerr <- read_bcerror(bcerr_path)
#' call_bcerror_sites(bcerr, min_error = 0.05, min_cov = 10)
call_bcerror_sites <- function(
  bcerror,
  min_error = 0.1,
  min_cov = 20,
  refs = NULL
) {
  required <- c("ref", "pos", "cov", "error_rate")
  missing <- setdiff(required, names(bcerror))
  if (length(missing) > 0) {
    cli_abort(
      "{.arg bcerror} is missing column{?s} {.field {missing}}."
    )
  }

  if (!is.null(refs)) {
    bcerror <- dplyr::filter(bcerror, .data$ref %in% refs)
  }

  bcerror |>
    dplyr::filter(
      .data$cov >= min_cov,
      .data$error_rate >= min_error
    ) |>
    dplyr::select("ref", "pos", "cov", "error_rate") |>
    dplyr::arrange(.data$ref, .data$pos)
}

#' Compute odds ratios between modification and charging.
#'
#' For each site, use individual reads as the unit of observation to ask
#' whether a read being modified at that site is associated with that same
#' read being charged. Each site yields a 2x2 table of modified/unmodified
#' against charged/uncharged, an odds ratio, and Fisher's exact test.
#'
#' An odds ratio above 1 means modified reads are more likely to be charged
#' than unmodified reads of the same tRNA; below 1 means the opposite. Because
#' both variables are measured on the same read, this is a within-tRNA
#' comparison and is not confounded by differences in abundance or charging
#' between tRNAs or between samples.
#'
#' Sites can come from either of two sources, set by `sites`:
#'
#' - Direct modification calls (`sites = NULL`, the default): every position
#'   present in `calls` is tested.
#' - Base-calling error: pass a site list from [call_bcerror_sites()] to
#'   restrict testing to positions called from error rates. In this case
#'   `calls` should be a per-read mismatch table rather than modification
#'   calls, so that the per-read status and the site selection come from the
#'   same signal.
#'
#' Reads with no call at a site are excluded from that site's table rather
#' than counted as unmodified, so `total_obs` varies between sites.
#'
#' @param calls Per-read calls, either a path to a `mod_calls.tsv.gz` /
#'   `mismatch_calls.tsv.gz` file or a tibble. Requires columns `read_id`,
#'   `chrom` (or `ref`), `ref_position` (or `pos`), and either `call_code`
#'   (where `"-"` means no call) or an integer `modified` column. A logical
#'   `within_alignment` column, if present, is used to drop calls outside the
#'   alignment.
#' @param charging Per-read charging status, either a path to a
#'   `charging_prob.tsv.gz` file or a tibble from [read_charging_calls()].
#' @param sites Optional tibble of sites to test, with columns `ref` and
#'   `pos`, typically from [call_bcerror_sites()]. If `NULL` (default), all
#'   positions in `calls` are tested.
#' @param refs Optional character vector of reference names to include. If
#'   `NULL` (default), all references are processed.
#' @param min_reads Minimum number of reads with both a call and a charging
#'   status required for a site to be tested. Default `10`.
#' @param ml_threshold Charging likelihood threshold, used only when
#'   `charging` is a path. Default `200`.
#' @param p_method Method for p-value adjustment, passed to [stats::p.adjust()].
#'   Default `"BH"`.
#'
#' @return A tibble with one row per tested site and columns: `ref`, `pos`,
#'   `n11` (modified and charged), `n10` (modified and uncharged), `n01`
#'   (unmodified and charged), `n00` (unmodified and uncharged), `total_obs`,
#'   `mod_freq`, `charged_freq`, `odds_ratio`, `log_odds_ratio`, `p_value`,
#'   and `p_adjusted`.
#'
#' @seealso [call_bcerror_sites()], [read_charging_calls()],
#'   [compute_odds_ratios()]
#'
#' @export
#'
#' @examples
#' calls <- clover_example("ecoli/mod_calls.tsv.gz")
#' charging <- clover_example("ecoli/charging_calls.tsv.gz")
#' compute_charging_odds_ratios(calls, charging, min_reads = 5)
#'
#' # Restrict to sites called from base-calling error instead.
#' bcerr <- read_bcerror(clover_example(
#'   "ecoli/summary/tables/wt-15-ctl-01/wt-15-ctl-01.bcerror.tsv.gz"
#' ))
#' sites <- call_bcerror_sites(bcerr, min_error = 0.05, min_cov = 10)
#' mismatches <- clover_example("ecoli/mismatch_calls.tsv.gz")
#' compute_charging_odds_ratios(
#'   mismatches,
#'   charging,
#'   sites = sites,
#'   min_reads = 5
#' )
compute_charging_odds_ratios <- function(
  calls,
  charging,
  sites = NULL,
  refs = NULL,
  min_reads = 10,
  ml_threshold = 200,
  p_method = "BH"
) {
  calls_tbl <- as_calls_table(calls)
  charging_tbl <- as_charging_table(charging, ml_threshold = ml_threshold)

  if (!is.null(refs)) {
    calls_tbl <- dplyr::filter(calls_tbl, .data$ref %in% refs)
  }

  if (!is.null(sites)) {
    calls_tbl <- restrict_to_sites(calls_tbl, sites)
  }

  # Charging is per read, and a read aligns to a single reference, so joining
  # on read_id alone is sufficient and avoids depending on the two files
  # naming that reference identically.
  calls_tbl <- dplyr::inner_join(
    calls_tbl,
    dplyr::select(charging_tbl, "read_id", "charged"),
    by = "read_id"
  )

  if (nrow(calls_tbl) == 0) {
    return(empty_charging_or())
  }

  results <- lapply(
    split(calls_tbl, calls_tbl$ref),
    charging_or_one_ref,
    min_reads = min_reads
  )

  out <- dplyr::bind_rows(results)

  if (nrow(out) == 0) {
    return(empty_charging_or())
  }

  out |>
    dplyr::mutate(
      p_adjusted = stats::p.adjust(.data$p_value, method = p_method)
    ) |>
    dplyr::arrange(.data$p_value, .data$ref, .data$pos)
}

# Compute site-versus-charging odds ratios for a single reference.
#
# `data` holds one row per (read, position) with `modified` and `charged`.
# Reads missing a call at a position are absent from the wide matrix (NA) and
# are excluded from that position's table rather than treated as unmodified.
charging_or_one_ref <- function(data, min_reads) {
  ref <- data$ref[1]

  wide <- data |>
    dplyr::summarize(
      modified = max(.data$modified),
      .by = c("read_id", "pos")
    ) |>
    tidyr::pivot_wider(
      names_from = "pos",
      values_from = "modified"
    )

  # Index rather than join: the charging vector has to line up with the rows of
  # the matrix built from `wide`, and a join does not guarantee that order.
  chg_map <- dplyr::distinct(data, .data$read_id, .data$charged)
  charged <- chg_map$charged[match(wide$read_id, chg_map$read_id)]

  pos_cols <- setdiff(names(wide), "read_id")
  if (length(pos_cols) == 0) {
    return(NULL)
  }

  mat <- as.matrix(wide[, pos_cols, drop = FALSE])

  is_mod <- !is.na(mat) & mat == 1
  is_unmod <- !is.na(mat) & mat == 0
  storage.mode(is_mod) <- "double"
  storage.mode(is_unmod) <- "double"

  chg <- as.numeric(charged)
  unchg <- 1 - chg

  n11 <- as.vector(crossprod(is_mod, chg))
  n10 <- as.vector(crossprod(is_mod, unchg))
  n01 <- as.vector(crossprod(is_unmod, chg))
  n00 <- as.vector(crossprod(is_unmod, unchg))

  total <- n11 + n10 + n01 + n00

  # Drop sites with too few observations, or with an empty margin: with no
  # variation in modification or in charging the odds ratio is undefined.
  keep <-
    total >= min_reads &
    (n11 + n10) > 0 &
    (n01 + n00) > 0 &
    (n11 + n01) > 0 &
    (n10 + n00) > 0

  if (!any(keep)) {
    return(NULL)
  }

  n11 <- n11[keep]
  n10 <- n10[keep]
  n01 <- n01[keep]
  n00 <- n00[keep]
  total <- total[keep]

  # Sample odds ratio, with a Haldane correction when any cell is empty. This
  # matches the definition used by compute_odds_ratios().
  any_zero <- n11 == 0 | n10 == 0 | n01 == 0 | n00 == 0
  odds_ratio <- ifelse(
    any_zero,
    ((n11 + 0.5) * (n00 + 0.5)) / ((n10 + 0.5) * (n01 + 0.5)),
    (n11 * n00) / (n10 * n01)
  )

  p_value <- vapply(
    seq_along(n11),
    function(i) {
      stats::fisher.test(
        matrix(c(n11[i], n01[i], n10[i], n00[i]), nrow = 2)
      )$p.value
    },
    numeric(1)
  )

  tibble::tibble(
    ref = ref,
    pos = as.integer(pos_cols[keep]),
    n11 = n11,
    n10 = n10,
    n01 = n01,
    n00 = n00,
    total_obs = total,
    mod_freq = (n11 + n10) / total,
    charged_freq = (n11 + n01) / total,
    odds_ratio = odds_ratio,
    log_odds_ratio = log(odds_ratio),
    p_value = p_value
  )
}

empty_charging_or <- function() {
  tibble::tibble(
    ref = character(),
    pos = integer(),
    n11 = numeric(),
    n10 = numeric(),
    n01 = numeric(),
    n00 = numeric(),
    total_obs = numeric(),
    mod_freq = numeric(),
    charged_freq = numeric(),
    odds_ratio = numeric(),
    log_odds_ratio = numeric(),
    p_value = numeric(),
    p_adjusted = numeric()
  )
}

# Normalize per-read calls (a path or a tibble, in modkit or mismatch layout)
# to columns read_id, ref, pos, modified.
as_calls_table <- function(calls) {
  if (is.character(calls)) {
    if (length(calls) != 1) {
      cli_abort("{.arg calls} must be a single path or a data frame.")
    }
    calls <- readr::read_tsv(calls, show_col_types = FALSE)
  }

  if (!is.data.frame(calls)) {
    cli_abort("{.arg calls} must be a single path or a data frame.")
  }

  if ("within_alignment" %in% names(calls)) {
    calls <- dplyr::filter(calls, .data$within_alignment == TRUE)
  }

  ref_col <- if ("chrom" %in% names(calls)) "chrom" else "ref"
  pos_col <- if ("ref_position" %in% names(calls)) "ref_position" else "pos"

  required <- c("read_id", ref_col, pos_col)
  missing <- setdiff(required, names(calls))
  if (length(missing) > 0) {
    cli_abort("{.arg calls} is missing column{?s} {.field {missing}}.")
  }

  if (!"modified" %in% names(calls) && !"call_code" %in% names(calls)) {
    cli_abort(
      "{.arg calls} needs a {.field call_code} or {.field modified} column."
    )
  }

  modified <- if ("modified" %in% names(calls)) {
    as.integer(calls$modified)
  } else {
    as.integer(calls$call_code != "-")
  }

  tibble::tibble(
    read_id = calls[["read_id"]],
    ref = calls[[ref_col]],
    pos = as.integer(calls[[pos_col]]),
    modified = modified
  )
}

# Normalize per-read charging (a path or a tibble) to read_id, ref, charged.
as_charging_table <- function(charging, ml_threshold = 200) {
  if (is.character(charging)) {
    if (length(charging) != 1) {
      cli_abort("{.arg charging} must be a single path or a data frame.")
    }
    return(read_charging_calls(charging, ml_threshold = ml_threshold))
  }

  if (!is.data.frame(charging)) {
    cli_abort("{.arg charging} must be a single path or a data frame.")
  }

  if (!"charged" %in% names(charging)) {
    if (!"charging_likelihood" %in% names(charging)) {
      cli_abort(
        "{.arg charging} needs a {.field charged} or
         {.field charging_likelihood} column."
      )
    }
    charging <- dplyr::mutate(
      charging,
      charged = as.integer(.data$charging_likelihood >= ml_threshold)
    )
  }

  if (!"read_id" %in% names(charging)) {
    cli_abort("{.arg charging} is missing column {.field read_id}.")
  }

  dplyr::select(charging, dplyr::any_of(c("read_id", "ref", "charged")))
}

# Restrict a calls table to a (ref, pos) site list. The calls table keeps all
# of its columns; `ref_col`/`pos_col` name its coordinate columns, which differ
# between the normalized layout and raw modkit output.
restrict_to_sites <- function(
  calls_tbl,
  sites,
  ref_col = "ref",
  pos_col = "pos"
) {
  if (!is.data.frame(sites)) {
    cli_abort("{.arg sites} must be a data frame with {.field ref} and
               {.field pos} columns.")
  }

  missing <- setdiff(c("ref", "pos"), names(sites))
  if (length(missing) > 0) {
    cli_abort("{.arg sites} is missing column{?s} {.field {missing}}.")
  }

  sites <- dplyr::distinct(
    dplyr::transmute(
      sites,
      ref = .data$ref,
      pos = as.integer(.data$pos)
    )
  )

  by <- stats::setNames(c("ref", "pos"), c(ref_col, pos_col))
  out <- dplyr::semi_join(calls_tbl, sites, by = by)

  if (nrow(out) == 0) {
    cli_warn(
      "No calls overlap the {nrow(sites)} supplied site{?s}; check that
       {.arg sites} and {.arg calls} use the same coordinate system."
    )
  }

  out
}
