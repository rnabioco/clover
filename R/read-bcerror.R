#' Read base-calling error (bcerror) TSV files.
#'
#' bcerror files are generated with [detectrms](https://github.com/rnabioco/detectrms-rs)
#'
#' @examples
#' bcerr_path <- clover_example("yeast/grande.bcerr.tsv.gz")
#' read_bcerror(bcerr_path)
#'
#' @param bcerr_path Path to bcerror file.
#'
#' @return A tibble
#' @export
read_bcerror <- function(bcerr_path) {

  col_names <- c(
    "ref", "pos", "base", "strand",
    "cov", "q_mean", "q_median", "q_std",
    "mis", "ins", "del", "ACGT"
  )

  readr::read_tsv(
    bcerr_path,
    col_names = col_names,
    skip = 1,
    show_col_types = FALSE
  ) |>
    dplyr::mutate(
      pos = as.integer(pos)
    ) |>
    tidyr::separate(
      ACGT,
      into = c("na", "nc", "ng", "nt"),
      sep = ":",
      convert = TRUE
    )

}
