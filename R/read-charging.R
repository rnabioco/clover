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
    stop("'paths' must be a named character vector.", call. = FALSE)
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
    stop("'paths' must be a named character vector.", call. = FALSE)
  }

  tbls <- lapply(names(paths), function(sid) {
    tbl <- read_odds_ratios(paths[[sid]])
    tbl$sample_id <- sid
    tbl
  })

  dplyr::bind_rows(tbls)
}
