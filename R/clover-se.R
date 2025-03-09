#' Create a SummarizedExperiment
#'
#' @param bam list of BAM files, names are sample names
#' @param pod5 list of POD5 files, names are samples names
#'
#' @param mods TSV file containing modification information
#' @param fa path to the FASTA reference used in alignment
#'
#' @examples
#' counts <- read_counts(clover_example("grande.counts.tsv.gz"))
#' bcerror <- read_bcerror(clover_example("grande.bcerror.tsv.gz"))
#' fasta <- read_fasta(clover_example("yeast/trna-ref.fa.gz"))
#'
#' CloverSE(counts, fasta, bcerror)
#'
#' @export
CloverSE <- function(counts, fasta, bcerror = NULL, ...) {
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
#' @export
read_mod_annotations <- function(mods) {
  readr::read_tsv(mods)
}
