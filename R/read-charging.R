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
  if (is.null(names(paths))) {
    cli_abort("{.arg paths} must be a named character vector.")
  }

  tbls <- lapply(names(paths), function(sid) {
    tbl <- read_charging(paths[[sid]])
    tbl$sample_id <- sid
    tbl
  })

  dplyr::bind_rows(tbls)
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
  if (is.null(names(paths))) {
    cli_abort("{.arg paths} must be a named character vector.")
  }

  tbls <- lapply(names(paths), function(sid) {
    tbl <- read_odds_ratios(paths[[sid]])
    tbl$sample_id <- sid
    tbl
  })

  dplyr::bind_rows(tbls)
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

    # Compute Fisher's test for each pair
    for (i in seq_along(pos_cols)[-length(pos_cols)]) {
      for (j in (i + 1):length(pos_cols)) {
        p1 <- pos_cols[i]
        p2 <- pos_cols[j]

        tbl <- table(
          factor(mat_data[[p1]], levels = c(0, 1)),
          factor(mat_data[[p2]], levels = c(0, 1))
        )

        if (any(rowSums(tbl) == 0) || any(colSums(tbl) == 0)) {
          next
        }

        ft <- tryCatch(stats::fisher.test(tbl), error = function(e) NULL)
        if (is.null(ft)) {
          next
        }

        results[[length(results) + 1]] <- tibble::tibble(
          ref = r,
          pos1 = p1,
          pos2 = p2,
          odds_ratio = ft$estimate,
          log_odds_ratio = log(ft$estimate + 1e-10),
          p_value = ft$p.value,
          total_obs = nrow(mat_data)
        )
      }
    }
  }

  dplyr::bind_rows(results)
}
