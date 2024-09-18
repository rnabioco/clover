#' Read base-calling error ("bcerror") TSV files.
#'
#' bcerror files are generated with [detectrms](https://github.com/rnabioco/detectrms-rs)
#'
#' `error_rate` is the sum of `mis`, `ins`, and `del`.
#'
#' @examples
#' bcerr_path <- clover_example("yeast/grande.bcerr.tsv.gz")
#' read_bcerror(bcerr_path)
#'
#' @param bcerr_path Path to bcerror file.
#' @param raw Logical, return raw data.
#'
#' @return A tibble
#' @export
read_bcerror <- function(bcerr_path, raw = FALSE) {
  col_names <- c(
    "ref", "pos", "base", "strand",
    "cov", "q_mean", "q_median", "q_std",
    "mis", "ins", "del", "ACGT"
  )

  raw_tbl <-
    read_tsv(
      bcerr_path,
      col_names = col_names,
      skip = 1,
      show_col_types = FALSE
    )

  if (raw) {
    return(raw_tbl)
  }

  raw_tbl |>
    mutate(
      pos = as.integer(pos),
      across(c(ref, base), as.factor),
      error_rate = mis + ins + del
    ) |>
    separate(
      ACGT,
      into = c("na", "nc", "ng", "nt"),
      sep = ":",
      convert = TRUE
    ) |>
    select(
      ref:q_std, error_rate, everything(),
      -strand
    )
}
