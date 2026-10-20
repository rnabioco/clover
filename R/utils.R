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

#' Normalise tRNA identifiers to one canonical form
#'
#' The readers, Sprinzl coordinate tables, [modomics_mods()] and
#' [plot_tRNA_structure()] spell the same tRNA differently
#' (`host-tRNA-Ala-GGC-1-1`, `nuc-tRNA-Ala-AGC-1-1`, `tRNA-Ile-GAU-1-1`,
#' `trna-Gln-TTG`). `normalize_trna_id()` maps them onto one form so that
#' tables can be joined on the result.
#'
#' The canonical form is `tRNA-<isotype>-<anticodon>`, with the anticodon in
#' DNA letters (T, not U), for example `tRNA-Ala-GGC`. It drops:
#'
#' - the source prefix: `host-`, `phage-`, `nuc-` or `mito-`;
#' - a trailing `-uncharged` (the reference variant in charging tables);
#' - the gene-copy suffix (`-1-1`), unless `copy = TRUE`.
#'
#' It also upper-cases the `tRNA` stem and converts U to T in the anticodon.
#' Isotype names such as `Ile2`, `fMet` and `SeC` are kept as written. This
#' is the same form [plot_tRNA_structure()] uses for its `trna` argument.
#'
#' Copy suffixes are not stable between sources: the same gene can be
#' `tRNA-Thr-CGT-1-2` in one table and `tRNA-Thr-CGT-2-1` in another. The
#' default therefore matches at the isodecoder level (isotype and anticodon).
#' Use `copy = TRUE` only when both tables come from the same source.
#'
#' The function is idempotent. A string that is not a tRNA identifier in
#' one of these forms (for example `"tRNA-Ala"`, which has no anticodon)
#' becomes `NA`, so `x[is.na(normalize_trna_id(x))]` lists what did not
#' resolve.
#'
#' @param x Character vector of tRNA identifiers.
#' @param copy If `TRUE`, keep the gene-copy suffix when `x` has one.
#'
#' @return Character vector the same length as `x`.
#'
#' @export
#' @examples
#' normalize_trna_id(c(
#'   "host-tRNA-Ala-GGC-1-1",
#'   "nuc-tRNA-Ala-AGC-1-1",
#'   "tRNA-Ile-GAU-1-1",
#'   "mito-tRNA-fMet-CAU",
#'   "tRNA-Gly-GCC-1-1-uncharged"
#' ))
#' normalize_trna_id("phage-tRNA-Thr-CGT-2-1", copy = TRUE)
normalize_trna_id <- function(x, copy = FALSE) {
  x <- as.character(x)
  if (length(x) == 0) {
    return(character())
  }
  parts <- stringr::str_match(
    x,
    paste0(
      "^(?:(?:host|phage|nuc|mito)-)?[tT][rR][nN][aA]-",
      "([A-Za-z0-9]+)-([ACGTUacgtu]{3})",
      "((?:-\\d+){0,2})(?:-uncharged)?$"
    )
  )
  suffix <- if (copy) parts[, 4] else ""
  out <- paste0(
    "tRNA-",
    parts[, 2],
    "-",
    chartr("U", "T", toupper(parts[, 3])),
    suffix
  )
  out[is.na(parts[, 1])] <- NA_character_
  out
}
