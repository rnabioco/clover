# Global tRNA coordinate functions ---------------------------------------------

#' Available organisms with pre-computed global coordinates
#'
#' @return Character vector of organism names
#' @export
#'
#' @examples
#' available_organisms()
available_organisms <- function() {
  c("ecoliK12", "sacCer", "hg38")
}

#' Load global tRNA coordinates
#'
#' Load pre-computed global coordinate mappings for tRNAs. These coordinates
#' enable alignment of tRNAs with different sequences onto a common axis
#' based on their structural positions (Sprinzl coordinates).
#'
#' @param organism Character string specifying the organism. One of
#'   "ecoliK12" (E. coli K12), "sacCer" (S. cerevisiae), or "hg38" (H. sapiens).
#' @param path Optional path to a custom coordinate file. If provided,
#'   `organism` is ignored.
#'
#' @return A tibble with columns:
#'   \describe{
#'     \item{trna_id}{tRNA identifier matching reference names}
#'     \item{seq_index}{Position in the tRNA sequence (1-based)}
#'     \item{sprinzl_index}{Sprinzl position number (1-76, or -1 for unlabeled)}
#'     \item{sprinzl_label}{Sprinzl label (e.g., "20", "20A", "e11")}
#'     \item{residue}{Nucleotide at this position}
#'     \item{global_index}{Global coordinate for cross-tRNA alignment (1..K)}
#'     \item{region}{Structural region (e.g., "D-loop", "anticodon-stem")}
#'   }
#'
#' @export
#'
#' @examples
#' coords <- load_global_coords("sacCer")
#' head(coords)
load_global_coords <- function(organism = "sacCer", path = NULL) {
  if (!is.null(path)) {
    return(read_coords_file(path))
  }

  organism <- match.arg(organism, available_organisms())


  coord_file <- system.file(
    "extdata", "coords",
    paste0(organism, "_global_coords.tsv.gz"),
    package = "clover",
    mustWork = TRUE
  )

  read_coords_file(coord_file)
}

#' Read coordinate file
#' @noRd
read_coords_file <- function(path) {
  readr::read_tsv(
    path,
    col_types = readr::cols(
      trna_id = readr::col_character(),
      source_file = readr::col_character(),
      seq_index = readr::col_integer(),
      sprinzl_index = readr::col_integer(),
      sprinzl_label = readr::col_character(),
      residue = readr::col_character(),
      sprinzl_ordinal = readr::col_double(),
      sprinzl_continuous = readr::col_double(),
      global_index = readr::col_integer(),
      region = readr::col_character()
    )
  ) |>
    dplyr::select(
      trna_id,
      seq_index,
      sprinzl_index,
      sprinzl_label,
      residue,
      global_index,
      region
    )
}

#' Add global coordinates to base-calling error data
#'
#' Join bcerror data with global coordinates to enable cross-tRNA comparisons.
#' Matches on tRNA ID and sequence position.
#'
#' @param bcerror Tibble of bcerror data from [read_bcerror()]
#' @param coords Tibble of global coordinates from [load_global_coords()],
#'   or a character string specifying an organism name.
#'
#' @return The input bcerror tibble with additional columns:
#'   `global_index`, `sprinzl_label`, and `region`.
#'
#' @export
#'
#' @examples
#' bcerr_path <- clover_example("yeast/grande.bcerr.tsv.gz")
#' bcerr <- read_bcerror(bcerr_path)
#' bcerr_with_coords <- add_global_coords(bcerr, "sacCer")
#' head(bcerr_with_coords)
add_global_coords <- function(bcerror, coords) {
  if (is.character(coords)) {
    coords <- load_global_coords(coords)
  }

  # Join on tRNA ID and position

  dplyr::left_join(
    bcerror,
    coords,
    by = c("ref" = "trna_id", "pos" = "seq_index")
  )
}

#' Get unique global index labels for plotting
#'
#' Returns a mapping of global indices to Sprinzl labels for axis annotation.
#'
#' @param coords Tibble of global coordinates from [load_global_coords()]
#'
#' @return A named character vector where names are global indices and
#'   values are Sprinzl labels.
#'
#' @export
#'
#' @examples
#' coords <- load_global_coords("sacCer")
#' labels <- get_global_labels(coords)
#' head(labels)
get_global_labels <- function(coords) {
  label_map <- coords |>
    dplyr::filter(sprinzl_label != "-1") |>
    dplyr::distinct(global_index, sprinzl_label) |>
    dplyr::arrange(global_index)

  stats::setNames(label_map$sprinzl_label, label_map$global_index)
}

#' Get structural regions for global indices
#'
#' Returns a mapping of global indices to structural regions for annotation.
#'
#' @param coords Tibble of global coordinates from [load_global_coords()]
#'
#' @return A tibble with columns `global_index`, `region`, `start`, and `end`
#'   defining contiguous regions.
#'
#' @export
#'
#' @examples
#' coords <- load_global_coords("sacCer")
#' regions <- get_region_bounds(coords)
#' regions
get_region_bounds <- function(coords) {
  coords |>
    dplyr::distinct(global_index, region) |>
    dplyr::arrange(global_index) |>
    dplyr::mutate(
      region_change = region != dplyr::lag(region, default = "")
    ) |>
    dplyr::mutate(
      region_group = cumsum(region_change)
    ) |>
    dplyr::group_by(region_group, region) |>
    dplyr::summarize(
      start = min(global_index),
      end = max(global_index),
      .groups = "drop"
    ) |>
    dplyr::select(-region_group)
}
