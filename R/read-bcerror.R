#' Read base-calling error (bcerror) TSV files.
#'
#' bcerror files can be created by x, y, z.
#'
#' @examples
#' bcerror_file <- clover_example("grande.bcerror.tsv.xz")
#' read_bcerror(bcerror_file)
#'
#' @return A tibble
#' @param path Path to bcerror file.
#' @export
read_bcerror <- function(path) {
  col_names <- c(
    "ref", "pos", "n_span", "n_base_mapped",
    "freq_a", "freq_t", "freq_g", "freq_c", "freq_n",
    "freq_mm", "freq_ins", "freq_del", "freq_bcerror",
    "qual_mean"
  )

  readr::read_tsv(path, col_names = col_names, skip = 1)
}
