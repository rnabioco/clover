#' Provide working directory for clover example files.
#'
#' @param path path to file
#'
#' @return A string with the path to the example file.
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

#' Shorten tRNA name strings for display.
#'
#' Strips source prefixes (e.g., `host-`, `phage-`), the `tRNA-` prefix,
#' and the gene-copy suffix (`-\d+-\d+`) to produce compact labels
#' suitable for plot axes.
#'
#' @param x Character vector of tRNA name strings.
#' @param strip_prefix Regex pattern for source prefixes to remove.
#'   Default `"^(host|phage)-"`.
#'
#' @return Character vector of shortened names (e.g.,
#'   `"host-tRNA-Ser-CGT-1-1"` becomes `"Ser-CGT"`).
#'
#' @export
#' @examples
#' shorten_trna_names("host-tRNA-Ser-CGT-1-1")
#' shorten_trna_names(c("phage-tRNA-Glu-TTC-2-1", "tRNA-Ala-AGC-3-1"))
shorten_trna_names <- function(x, strip_prefix = "^(host|phage)-") {
  x |>
    sub(strip_prefix, "", x = _) |>
    sub("^tRNA-", "", x = _) |>
    sub("-\\d+-\\d+$", "", x = _)
}
