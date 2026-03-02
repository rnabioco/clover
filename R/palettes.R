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
    "Ala" = "#D5F6E2",
    "Arg" = "#A5C2F6",
    "Asn" = "#FF9902",
    "Asp" = "#9A00FF",
    "Cys" = "#D9B5FF",
    "Gln" = "#44B3E1",
    "Glu" = "#D10C0D",
    "Gly" = "#000000",
    "His" = "#0070C0",
    "Ile" = "#FFFF00",
    "Leu" = "#FF00FF",
    "Lys" = "#FFD9AF",
    "Met" = "#ADADAD",
    "Phe" = "#C99D89",
    "Pro" = "#12501A",
    "Ser" = "#FFC000",
    "Thr" = "#002060",
    "Trp" = "#800000",
    "Tyr" = "#47D359",
    "Val" = "#FFCAE2",
    "fMet" = "#ADADAD",
    "Ile2" = "#FFFF00",
    "SeC" = "#FFC000"
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

#' @noRd
nucleotide_colors <- function() {
  c(
    "A" = "#4DAF4A",
    "C" = "#377EB8",
    "G" = "#FFD92F",
    "U" = "#E41A1C"
  )
}

#' @noRd
region_colors <- function() {
  c(
    "acceptor-stem" = "#E41A1C",
    "acceptor-tail" = "#E41A1C",
    "D-stem" = "#377EB8",
    "D-loop" = "#4DAF4A",
    "anticodon-stem" = "#984EA3",
    "anticodon-loop" = "#FF7F00",
    "variable-region" = "#A65628",
    "variable-arm" = "#A65628",
    "T-stem" = "#F781BF",
    "T-loop" = "#999999",
    "unknown" = "grey70"
  )
}
