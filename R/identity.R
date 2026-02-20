# Identity element functions ---------------------------------------------------

#' Retrieve tRNA aminoacylation identity elements
#'
#' Returns experimentally validated identity elements (determinants
#' and antideterminants) for tRNA aminoacylation, based on data from
#' Giege & Eriani (2023). Elements are defined per amino acid family
#' using Sprinzl numbering.
#'
#' @param organism Character string specifying the organism name
#'   (e.g., `"Escherichia coli"`, `"Saccharomyces cerevisiae"`,
#'   `"Homo sapiens"`). Use [identity_organisms()] to list supported
#'   organisms.
#' @param amino_acid Optional character vector of 3-letter amino acid
#'   codes to filter by (e.g., `"Ala"`, `"Arg"`). If `NULL` (default),
#'   returns elements for all amino acids.
#' @param type Which element types to return: `"both"` (default),
#'   `"determinant"`, or `"antideterminant"`.
#'
#' @return A tibble with columns:
#'   - `amino_acid`: 3-letter amino acid code
#'   - `domain`: evolutionary domain (`"Bacteria"` or `"Eukarya"`)
#'   - `sprinzl_pos`: Sprinzl position number (integer; `NA` for
#'     structural determinants like the long variable arm; `-1` for
#'     the His G\eqn{_{-1}} position)
#'   - `nucleotide`: expected base (A/C/G/U, or `NA` if variable)
#'   - `region`: tRNA structural region name
#'   - `type`: `"determinant"` or `"antideterminant"`
#'   - `strength`: `"strong"` or `"weak"` (determinants only)
#'   - `pair_pos`: partner Sprinzl position for base pairs (`NA`
#'     for single-nucleotide elements)
#'   - `pair_type`: `"wc"` or `"wobble"` (`NA` for unpaired)
#'   - `against_aars`: aaRS blocked (antideterminants only)
#'   - `universal`: `TRUE` if conserved across all domains
#'   - `description`: human-readable description
#'
#' @export
#'
#' @references
#' Giege R, Eriani G (2023). "The tRNA identity landscape for
#' aminoacylation and beyond." *Nucleic Acids Research*, 51(4),
#' 1528--1570. \doi{10.1093/nar/gkad007}
#'
#' @examples
#' # All identity elements for E. coli
#' identity_elements("Escherichia coli")
#'
#' # Ala determinants only
#' identity_elements("Escherichia coli", amino_acid = "Ala",
#'   type = "determinant")
identity_elements <- function(
  organism,
  amino_acid = NULL,
  type = c("both", "determinant", "antideterminant")
) {
  type <- match.arg(type)

  org_map <- load_identity_organism_map()
  domain <- org_map$domain[org_map$organism == organism]

  if (length(domain) == 0) {
    supported <- paste0(
      "'",
      org_map$organism,
      "'",
      collapse = ", "
    )
    cli::cli_abort(
      c(
        "Organism {.val {organism}} is not supported.",
        "i" = "Supported organisms: {supported}.",
        "i" = "Use {.fn identity_organisms} to list them."
      )
    )
  }

  result <- dplyr::tibble(
    amino_acid = character(),
    domain = character(),
    sprinzl_pos = integer(),
    nucleotide = character(),
    region = character(),
    type = character(),
    strength = character(),
    pair_pos = integer(),
    pair_type = character(),
    against_aars = character(),
    universal = logical(),
    description = character()
  )

  if (type %in% c("both", "determinant")) {
    dets <- load_cached_determinants()
    dets <- dets[dets$domain == domain, ]
    dets$type <- "determinant"
    dets$against_aars <- NA_character_
    result <- dplyr::bind_rows(result, dets)
  }

  if (type %in% c("both", "antideterminant")) {
    antis <- load_cached_antideterminants()
    antis <- antis[antis$domain == domain, ]
    antis$type <- "antideterminant"
    antis$strength <- NA_character_
    antis$universal <- NA
    result <- dplyr::bind_rows(result, antis)
  }

  result <- result |>
    dplyr::select(
      amino_acid,
      domain,
      sprinzl_pos,
      nucleotide,
      region,
      type,
      strength,
      pair_pos,
      pair_type,
      against_aars,
      universal,
      description
    )

  if (!is.null(amino_acid)) {
    result <- result[result$amino_acid %in% amino_acid, ]
  }

  result
}

#' List supported organisms for identity elements
#'
#' Returns the names of organisms for which tRNA aminoacylation
#' identity element data is available.
#'
#' @return A character vector of organism names.
#'
#' @export
#'
#' @examples
#' identity_organisms()
identity_organisms <- function() {
  org_map <- load_identity_organism_map()
  org_map$organism
}

#' Map identity elements to tRNA sequence positions
#'
#' Converts Sprinzl positions in identity element data to 1-based
#' sequence positions for a specific tRNA, enabling overlay on
#' structure plots via the `outlines` parameter of
#' [plot_tRNA_structure()].
#'
#' @param elements A tibble of identity elements as returned by
#'   [identity_elements()].
#' @param sprinzl_coords A tibble of Sprinzl coordinates as returned
#'   by [read_sprinzl_coords()].
#' @param trna_id Character string identifying the tRNA to map
#'   (must match a `trna_id` value in `sprinzl_coords`).
#'
#' @return The input `elements` tibble with an added `pos` column
#'   containing the 1-based sequence position. Rows where the
#'   Sprinzl position does not exist in the target tRNA (including
#'   His G\eqn{_{-1}} at `sprinzl_pos = -1` and structural
#'   determinants with `sprinzl_pos = NA`) are dropped.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' coords <- read_sprinzl_coords(
#'   clover_example("sprinzl/sacCer_global_coords.tsv.gz")
#' )
#' elems <- identity_elements("Saccharomyces cerevisiae",
#'   amino_acid = "Ala")
#' map_identity_to_trna(elems, coords, "nuc-tRNA-Ala-AGC-1-1")
#' }
map_identity_to_trna <- function(elements, sprinzl_coords, trna_id) {
  coords <- sprinzl_coords[sprinzl_coords$trna_id == trna_id, ]

  if (nrow(coords) == 0) {
    cli::cli_abort(
      "tRNA {.val {trna_id}} not found in {.arg sprinzl_coords}."
    )
  }

  # Convert sprinzl_pos to character for joining with sprinzl_label
  elements$sprinzl_int <- as.character(elements$sprinzl_pos)

  merged <- dplyr::inner_join(
    elements,
    coords[, c("sprinzl_label", "pos")],
    by = c("sprinzl_int" = "sprinzl_label")
  )

  merged$sprinzl_int <- NULL
  merged
}

# Internal helpers -------------------------------------------------------------

load_cached_determinants <- function() {
  path <- system.file(
    "extdata",
    "identity",
    "determinants.rds",
    package = "clover"
  )
  readRDS(path)
}

load_cached_antideterminants <- function() {
  path <- system.file(
    "extdata",
    "identity",
    "antideterminants.rds",
    package = "clover"
  )
  readRDS(path)
}

load_identity_organism_map <- function() {
  path <- system.file(
    "extdata",
    "identity",
    "organism_map.rds",
    package = "clover"
  )
  readRDS(path)
}
