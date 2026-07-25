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
#' # Pruning sites before testing
#'
#' A genome-wide site list is mostly sites that were never testable, and each
#' one still consumes FDR budget. Two arguments drop them up front. `min_margin`
#' requires a minimum count of modified, unmodified, charged and uncharged
#' reads; `max_p` drops sites whose margins make the intended significance
#' threshold unreachable no matter how the reads fall.
#'
#' Both filter on the margins of the 2x2 table, never on the observed odds
#' ratio, and that distinction matters. The margins are ancillary, so filtering
#' on them leaves the null distribution of the surviving p-values intact.
#' Filtering on effect size would enrich for small p-values and invalidate the
#' correction applied afterwards, so do not pre-filter `sites` on an effect
#' measured from the same reads.
#'
#' # Sites near the 3' end are not interpretable
#'
#' The amino acid is esterified to the terminal A, and charging is called from
#' the nanopore signal over that CCA end. A site there is therefore not
#' independent of the charging call it is being tested against: the base caller
#' is reading out the same thing the charging model is, and the odds ratio
#' between them is close to tautological.
#'
#' This dominates the result. Across two zebrafish samples, 82% of terminal-A
#' sites came back significant with a median log odds ratio of 1.9, and the
#' effect decayed monotonically with distance from that residue -- roughly 50%
#' of sites one or two nucleotides in, 16-33% from three to ten, and 4.6%
#' beyond fifty, where it flattens out. By region, `cca` and `discriminator`
#' ran at 63% and 36% against 3-5% in the D-loop and anticodon arms.
#'
#' Note that the position number of the terminal A varies between references,
#' since tRNAs differ in length, so filter on distance from the 3' end or on
#' the `region` column of a Sprinzl coordinate table rather than on an absolute
#' position. Dropping `cca`, `discriminator` and the 3' acceptor stem removes
#' the worst of it; the gradient suggests discarding everything within about
#' twenty nucleotides of the end. This caveat applies whichever site source is
#' used, since it concerns the charging call rather than the modification call.
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
#' @param min_margin Minimum count in each margin of the 2x2 table: modified,
#'   unmodified, charged and uncharged reads. Default `1`, which only requires
#'   the odds ratio to be defined. Raising it drops sites too thin to carry
#'   power, such as a position where three of ten thousand reads are modified.
#' @param max_p Skip sites where no arrangement of the reads could reach this
#'   p-value, given the margins. Default `1`, which tests everything. Set it to
#'   the significance threshold you intend to use.
#' @param ml_threshold Charging likelihood threshold, used only when
#'   `charging` is a path. Default `200`.
#' @param dedupe Whether to collapse repeated calls for the same read and
#'   position before counting. Needed for `modkit` output, where a read with
#'   several modification channels at one position contributes a row per
#'   channel. Default `TRUE`. Set `FALSE` for per-read mismatch calls, which
#'   carry one row per read and position, to skip a grouping pass that is
#'   costly on large tables.
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
  min_margin = 1,
  max_p = 1,
  ml_threshold = 200,
  dedupe = TRUE,
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

  if (dedupe) {
    calls_tbl <- dplyr::summarize(
      calls_tbl,
      modified = max(.data$modified),
      .by = c("read_id", "ref", "pos", "charged")
    )
  }

  ref_f <- factor(calls_tbl$ref)
  pos_f <- factor(calls_tbl$pos, levels = sort(unique(calls_tbl$pos)))

  out <- charging_odds_ratios_cpp(
    as.integer(ref_f),
    as.integer(pos_f),
    as.integer(calls_tbl$modified),
    as.integer(calls_tbl$charged),
    nlevels(ref_f),
    nlevels(pos_f),
    as.integer(min_reads),
    as.integer(min_margin),
    as.numeric(max_p)
  )

  if (nrow(out) == 0) {
    return(empty_charging_or())
  }

  tibble::tibble(
    ref = levels(ref_f)[out$ref_idx],
    pos = as.integer(levels(pos_f)[out$pos_idx]),
    n11 = out$n11,
    n10 = out$n10,
    n01 = out$n01,
    n00 = out$n00,
    total_obs = out$total_obs,
    mod_freq = (out$n11 + out$n10) / out$total_obs,
    charged_freq = (out$n11 + out$n01) / out$total_obs,
    odds_ratio = out$odds_ratio,
    log_odds_ratio = out$log_odds_ratio,
    p_value = out$p_value
  ) |>
    dplyr::mutate(
      p_adjusted = stats::p.adjust(.data$p_value, method = p_method)
    ) |>
    dplyr::arrange(.data$p_value, .data$ref, .data$pos)
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
