# Color palette functions ------------------------------------------------------

#' Amino acid color palette.
#'
#' Return a named character vector of colors for the 20 standard amino
#' acids plus fMet, Ile2, and SeC.
#'
#' @return A named character vector of hex colors.
#'
#' @export
#'
#' @examples
#' aa_colors()
aa_colors <- function() {
  c(
    "Ala" = "#1f77b4",
    "Arg" = "#ff7f0e",
    "Asn" = "#2ca02c",
    "Asp" = "#d62728",
    "Cys" = "#9467bd",
    "Gln" = "#8c564b",
    "Glu" = "#e377c2",
    "Gly" = "#7f7f7f",
    "His" = "#bcbd22",
    "Ile" = "#17becf",
    "Leu" = "#aec7e8",
    "Lys" = "#ffbb78",
    "Met" = "#98df8a",
    "Phe" = "#ff9896",
    "Pro" = "#c5b0d5",
    "Ser" = "#c49c94",
    "Thr" = "#f7b6d2",
    "Trp" = "#c7c7c7",
    "Tyr" = "#dbdb8d",
    "Val" = "#9edae5",
    "fMet" = "#393b79",
    "Ile2" = "#5254a3",
    "SeC" = "#6b6ecf"
  )
}

#' Charging status color palette.
#'
#' Return a named character vector of colors for charged and uncharged
#' tRNA states.
#'
#' @return A named character vector of hex colors.
#'
#' @export
#'
#' @examples
#' charging_colors()
charging_colors <- function() {
  c(
    "Charged" = "grey10",
    "Uncharged" = "#da7400ff"
  )
}
