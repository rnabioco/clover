#' Read a bedMethyl file from modkit pileup
#'
#' Read per-position modification percentages from a bedMethyl file
#' produced by `modkit pileup`. The 18-column bedMethyl format includes
#' per-position counts and percentages for each modification type.
#'
#' @param path Path to a bedMethyl file (may be gzipped).
#' @param min_cov Minimum valid coverage to retain a position. Default `1`.
#'
#' @return A tibble with columns: `ref`, `pos`, `strand`, `mod_code`,
#'   `n_valid`, `percent_mod`, `n_mod`, `n_canonical`.
#'
#' @export
#'
#' @examples
#' # read_bedmethyl("sample.bedmethyl.gz")
read_bedmethyl <- function(path, min_cov = 1) {
  bedmethyl_cols <- c(
    "chrom", "start", "end", "name", "score", "strand",
    "thickStart", "thickEnd", "itemRgb", "n_valid",
    "percent_mod", "n_mod", "n_canonical", "n_other_mod",
    "n_delete", "n_fail", "n_diff", "n_nocall"
  )

  raw <- readr::read_tsv(
    path,
    col_names = bedmethyl_cols,
    col_types = readr::cols(
      chrom = readr::col_character(),
      start = readr::col_integer(),
      end = readr::col_integer(),
      name = readr::col_character(),
      score = readr::col_integer(),
      strand = readr::col_character(),
      thickStart = readr::col_integer(),
      thickEnd = readr::col_integer(),
      itemRgb = readr::col_character(),
      n_valid = readr::col_integer(),
      percent_mod = readr::col_double(),
      n_mod = readr::col_integer(),
      n_canonical = readr::col_integer(),
      n_other_mod = readr::col_integer(),
      n_delete = readr::col_integer(),
      n_fail = readr::col_integer(),
      n_diff = readr::col_integer(),
      n_nocall = readr::col_integer()
    ),
    comment = "#",
    show_col_types = FALSE
  )

  raw |>
    dplyr::filter(.data$n_valid >= min_cov) |>
    dplyr::transmute(
      ref = .data$chrom,
      pos = .data$start,
      strand = .data$strand,
      mod_code = .data$name,
      n_valid = .data$n_valid,
      percent_mod = .data$percent_mod,
      n_mod = .data$n_mod,
      n_canonical = .data$n_canonical
    )
}

#' Summarize per-position modification frequency from modkit extract
#'
#' Read a `mod_calls.tsv.gz` file from `modkit extract` and compute
#' per-position modification frequency (proportion of reads with a
#' non-canonical call code).
#'
#' @param path Path to a `mod_calls.tsv.gz` file.
#' @param refs Optional character vector of reference names to include.
#' @param min_reads Minimum number of reads at a position to retain.
#'   Default `1`.
#'
#' @return A tibble with columns: `ref`, `pos`, `n_reads`, `n_mod`,
#'   `mod_freq`.
#'
#' @export
#'
#' @examples
#' path <- clover_example("ecoli/mod_calls.tsv.gz")
#' summarize_mod_calls(path)
summarize_mod_calls <- function(path, refs = NULL, min_reads = 1) {
  mc <- readr::read_tsv(path, show_col_types = FALSE) |>
    dplyr::filter(.data$within_alignment == TRUE)

  if (!is.null(refs)) {
    mc <- mc |> dplyr::filter(.data$chrom %in% refs)
  }

  mc |>
    dplyr::mutate(modified = as.integer(.data$call_code != "-")) |>
    dplyr::summarize(
      n_reads = dplyr::n(),
      n_mod = sum(.data$modified),
      .by = c("chrom", "ref_position")
    ) |>
    dplyr::filter(.data$n_reads >= min_reads) |>
    dplyr::transmute(
      ref = .data$chrom,
      pos = .data$ref_position,
      n_reads = .data$n_reads,
      n_mod = .data$n_mod,
      mod_freq = .data$n_mod / .data$n_reads
    )
}
