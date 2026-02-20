# Sprinzl coordinate functions ------------------------------------------------

#' Read Sprinzl coordinates from a global coordinates TSV file.
#'
#' Read a tRNAs-in-space global coordinates file and return a tibble with
#' standardized column names. These files map each position of each tRNA
#' to its Sprinzl numbering, structural region, and global alignment index.
#'
#' @param path Path to a global coordinates TSV file (may be gzipped).
#'
#' @return A tibble with columns:
#'   - `trna_id`: tRNA identifier (e.g., `nuc-tRNA-Ala-AGC-1-1`)
#'   - `pos`: 1-based position in tRNA body
#'   - `sprinzl_label`: Sprinzl position (character: "1", "20a", "47:e1", etc.)
#'   - `global_index`: cross-tRNA universal alignment position
#'   - `region`: structural region (acceptor-stem, D-stem, D-loop, etc.)
#'   - `residue`: reference nucleotide
#'
#' @export
#'
#' @examples
#' path <- clover_example("sprinzl/sacCer_global_coords.tsv.gz")
#' coords <- read_sprinzl_coords(path)
#' coords
read_sprinzl_coords <- function(path) {
  readr::read_tsv(
    path,
    col_types = readr::cols(sprinzl_label = readr::col_character()),
    show_col_types = FALSE
  ) |>
    dplyr::rename(pos = seq_index) |>
    dplyr::select(
      trna_id,
      pos,
      sprinzl_label,
      global_index,
      region,
      residue
    )
}

#' Order Sprinzl position labels.
#'
#' Create an ordered factor from Sprinzl position labels. Handles numeric
#' positions and lettered insertions (e.g., "20a", "20A", "47:e1").
#'
#' @param labels Character vector of Sprinzl labels.
#'
#' @return A factor with levels ordered by position number then suffix.
#'
#' @export
#'
#' @examples
#' order_sprinzl_positions(c("20a", "20", "1", "21"))
order_sprinzl_positions <- function(labels) {
  ordering <- tibble::tibble(label = unique(labels)) |>
    dplyr::filter(!is.na(label), label != "") |>
    dplyr::mutate(
      base_num = as.numeric(stringr::str_extract(label, "^\\d+")),
      suffix = stringr::str_extract(label, "[a-zA-Z]$|:e\\d+$"),
      suffix_order = dplyr::case_when(
        is.na(suffix) ~ 0,
        stringr::str_detect(suffix, "^[a-z]$") ~ match(suffix, letters),
        stringr::str_detect(suffix, "^[A-Z]$") ~ match(suffix, LETTERS),
        stringr::str_detect(suffix, "^:e") ~
          100 + as.numeric(stringr::str_extract(suffix, "\\d+$")),
        .default = 0
      ),
      sort_key = base_num + suffix_order / 1000
    ) |>
    dplyr::arrange(sort_key)

  factor(labels, levels = ordering$label)
}

#' tRNA structural regions.
#'
#' Return a named list mapping canonical tRNA structural region names
#' to integer vectors of Sprinzl position numbers.
#'
#' @return A named list of integer vectors.
#'
#' @export
#'
#' @examples
#' trna_regions()
trna_regions <- function() {
  list(
    acceptor_stem = c(1L:7L, 66L:72L),
    d_arm = 8L:25L,
    anticodon_stem = c(26L:31L, 39L:44L),
    anticodon_loop = 32L:38L,
    variable_loop = 44L:48L,
    t_arm = 49L:65L,
    discriminator = 73L,
    cca = 74L:76L
  )
}
