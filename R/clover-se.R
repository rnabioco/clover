#' Create a CloverExperiment object
#'
#' Creates a SummarizedExperiment-based object for storing tRNA sequencing data,
#' including expression counts, base-calling error rates, and global coordinates.
#'
#' @param counts Matrix or data.frame of tRNA counts (rows = tRNAs, cols = samples)
#' @param bcerror Tibble of base-calling error data from [read_bcerror()]
#' @param coords Tibble of global coordinates from [load_global_coords()],
#'   or a character string specifying an organism name.
#' @param metadata Named list of additional metadata (e.g., file paths)
#'
#' @return A CloverExperiment object
#'
#' @export
CloverSE <- function(
    counts = NULL,
    bcerror = NULL,
    coords = NULL,
    metadata = list()
) {
  # Handle coordinates
  if (is.character(coords)) {
    coords <- load_global_coords(coords)
  }

 # Add coordinates to bcerror if both provided
  if (!is.null(bcerror) && !is.null(coords)) {
    bcerror <- add_global_coords(bcerror, coords)
  }

  # Store in metadata
  metadata$bcerror <- bcerror
  metadata$coords <- coords

  # Create minimal SummarizedExperiment
  if (!is.null(counts)) {
    se <- SummarizedExperiment::SummarizedExperiment(
      assays = list(counts = as.matrix(counts)),
      metadata = metadata
    )
  } else {
    se <- SummarizedExperiment::SummarizedExperiment(
      metadata = metadata
    )
  }

  .CloverSE(se)
}

.CloverSE <- setClass("CloverExperiment", contains = "SummarizedExperiment")

#' Read FASTA reference
#'
#' @examples
#' fa <- clover_example("yeast/trna-ref.fa.gz")
#' read_fasta(fa)
#'
#' @param fa path to fasta file
#'
#' @import Biostrings
#' @export
read_fasta <- function(fa) {
  Biostrings::readDNAStringSet(fa)
}

#' Read modifications file
#'
#' TSV file with 4 columns:
#'
#' 1. `ref`
#' 2. `pos`
#' 3. `mod_full`
#' 4. `mod1`
#'
#' @param mods path to modifications file
#' @export
read_mod_annotations <- function(mods) {
  readr::read_tsv(mods)
}
