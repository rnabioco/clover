#' Create a SummarizedExperiment from pipeline output.
#'
#' Read a pipeline configuration file and load all result types into a
#' [SummarizedExperiment::SummarizedExperiment] object. The returned object
#' contains count matrices as assays, tRNA metadata as rowData, and
#' sample metadata as colData.
#'
#' @param config_path Path to a pipeline `config.yaml` file.
#' @param types Character vector of result types to load. Valid values:
#'   `"charging"`, `"bcerror"`, `"odds_ratios"`.
#' @param sample_info An optional data frame with a `sample_id` column and
#'   additional experimental factor columns (e.g., `condition`, `replicate`).
#'   If `NULL`, a minimal colData is created from sample names.
#' @param min_count Minimum total count across all samples for a tRNA
#'   to be retained in count matrices. Default `10`.
#'
#' @return A [SummarizedExperiment::SummarizedExperiment] with:
#'   - **assay "counts"**: abundance count matrix (charged + uncharged)
#'   - **assay "charging"**: charging count matrix (charged/uncharged columns
#'     per sample), only if `"charging"` is in `types`
#'   - **colData**: sample metadata
#'   - **metadata**: list with `$config`, `$bcerror`, `$odds_ratios`, `$fasta`
#'     as available
#'
#' @export
#'
#' @examples
#' \dontrun{
#' se <- create_clover("path/to/config.yaml")
#' SummarizedExperiment::assay(se, "counts")
#' SummarizedExperiment::colData(se)
#' }
create_clover <- function(
    config_path,
    types = c("charging", "bcerror", "odds_ratios"),
    sample_info = NULL,
    min_count = 10) {
  types <- match.arg(types, several.ok = TRUE)

  config <- read_pipeline_config(config_path)
  results <- read_pipeline_results(config_path, types = types)

  # Build assays from charging data
  assay_list <- list()

  if ("charging" %in% types && !is.null(results$charging)) {
    assay_list[["counts"]] <- abundance_count_matrix(
      results$charging,
      min_count = min_count
    )
  }

  if (length(assay_list) == 0) {
    stop("No assay data could be loaded. Check that charging data exists.",
      call. = FALSE
    )
  }

  # Build colData
  count_mat <- assay_list[["counts"]]
  coldata <- build_coldata(count_mat, sample_info)

  # Read reference FASTA if available
  fasta <- NULL
  if (!is.null(config$fasta) && file.exists(config$fasta)) {
    fasta <- tryCatch(
      read_fasta(config$fasta),
      error = function(e) NULL
    )
  }

  # Build rowData from tRNA names
  trna_names <- rownames(count_mat)
  row_data <- S4Vectors::DataFrame(
    tRNA = trna_names,
    row.names = trna_names
  )

  # Add sequence lengths from FASTA if available

  if (!is.null(fasta)) {
    fasta_names <- names(fasta)
    matched <- match(trna_names, fasta_names)
    row_data$seq_length <- ifelse(
      is.na(matched),
      NA_integer_,
      Biostrings::width(fasta)[matched]
    )
  }

  # Build metadata list
  meta <- list(config = config)
  if ("bcerror" %in% types && !is.null(results$bcerror)) {
    meta$bcerror <- results$bcerror
  }
  if ("odds_ratios" %in% types && !is.null(results$odds_ratios)) {
    meta$odds_ratios <- results$odds_ratios
  }
  if (!is.null(fasta)) {
    meta$fasta <- fasta
  }

  se <- SummarizedExperiment::SummarizedExperiment(
    assays = assay_list,
    rowData = row_data,
    colData = coldata,
    metadata = meta
  )

  se
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
#' @return A [Biostrings::DNAStringSet].
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
#' @return A tibble.
#' @export
read_mod_annotations <- function(mods) {
  readr::read_tsv(mods)
}
