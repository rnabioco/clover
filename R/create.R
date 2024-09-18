#' Create a RangedSummarizedExperiment object
#'
#' @param bam list of BAM files, names are sample names
#' @param pod5 list of POD5 files, names are samples names
#'
#' @param mods TSV file containing modification information
#' @param fa path to the FASTA reference used in alignment
#'
#' @importFrom SummarizedExperiment SummarizedExperiment
#' @importFrom GenomicRanges GRangesList
#'
#' @export
create_clover <- function(fa, bam, mods = NULL, pod5 = NULL) {
  SummarizedExperiment(
    assays = list(
      counts = list(),
      bcerror = list(),
    ),
    rowData = mods,
    rowRanges = GRangesList(),
    metadata = list(
      pod5 = pod5
    )
  )
}

#' @importFrom Biostrings readDNAStringSet
read_fasta <- function(fa) {
  ss <- readDNAStringSet(fa)
}

read_mods <- function(mods) {
  readr::read_tsv(mods)
}

# read_pod5 <- function(pod5) {
#  rpod5::read_pod5(pod5)
# }
