#' Provide working directory for clover example files.
#'
#' @param path path to file
#'
#' @examples
#' clover_example("ecoli/config.yaml")
#'
#' @export
clover_example <- function(path) {
  system.file("extdata", path, package = "clover", mustWork = TRUE)
}

#' Convert DNA anticodon to RNA in tRNA name strings
#'
#' Replaces T with U in the anticodon portion of tRNA identifiers,
#' e.g., `"tRNA-Glu-TTC-1-1"` becomes `"tRNA-Glu-UUC-1-1"`.
#'
#' @param x Character vector of tRNA name strings.
#'
#' @return Character vector with DNA anticodons converted to RNA.
#'
#' @export
#' @examples
#' dna_to_rna_anticodon("tRNA-Glu-TTC-1-1")
dna_to_rna_anticodon <- function(x) {
  parts <- stringr::str_match(x, "^(tRNA-[A-Za-z0-9]+-)(\\w+)(-.*)")
  ifelse(
    is.na(parts[, 1]),
    x,
    paste0(parts[, 2], chartr("T", "U", parts[, 3]), parts[, 4])
  )
}
