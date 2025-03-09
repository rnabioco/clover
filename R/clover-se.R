#' Create a SummarizedExperiment
#'
#' @param counts path to counts file
#' @param fasta path to the FASTA reference used in alignment
#' @param bcerror path to base-calling error file
#'
#' @param bam list of BAM files, names are sample names
#' @param pod5 list of POD5 files, names are samples names
#' @param mods TSV file containing modification information
#'
#' @export
CloverSE <- function(
  counts,
  fasta,
  bcerror = NULL,
  bam = NULL,
  pod5 = NULL,
  mods = NULL
) {
  se <- SummarizedExperiment::SummarizedExperiment(
    list(
      counts = counts,
      fasta = fasta,
      berror = bcerror,
      ...
    )
  )

  #  assays = list(
  #    counts = list(),
  #    bcerror = list()
  #  ),
  #  rowData = list(),
  #  rowRanges = GRangesList(),
  #  metadata = list(
  #    pod5 = list()
  #  )
  # )

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
