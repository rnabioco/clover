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

#' List available coordinate groups for an organism
#'
#' Each organism has multiple coordinate groups based on tRNA structural type
#' (Type I vs Type II) and labeling offset. This function lists all available
#' groups with their tRNA counts.
#'
#' @param organism Character string specifying the organism.
#'
#' @return A tibble with columns: `organism`, `offset`, `type`, `n_trnas`, `file`
#'
#' @export
#'
#' @examples
#' available_coord_groups("sacCer")
available_coord_groups <- function(organism = "sacCer") {
  organism <- match.arg(organism, available_organisms())

  # Find all coordinate files for this organism
  coord_dir <- system.file("extdata", "coords", package = "clover", mustWork = TRUE)
  pattern <- paste0("^", organism, "_global_coords_offset[+-]?[0-9]+_type[12]\\.tsv\\.gz$")
  files <- list.files(coord_dir, pattern = pattern, full.names = TRUE)

  if (length(files) == 0) {
    stop("No coordinate files found for organism: ", organism)
  }

  # Parse filenames and count tRNAs
  purrr::map_dfr(files, function(f) {
    meta <- parse_coord_filename(f)
    # Quick count of unique tRNAs without loading full file
    n_trnas <- length(unique(readr::read_tsv(f, col_types = "c------", col_select = 1)[[1]]))
    tibble::tibble(
      organism = meta$organism,
      offset = meta$offset,
      type = meta$type,
      n_trnas = n_trnas,
      file = basename(f)
    )
  }) |>
    dplyr::arrange(type, offset)
}

#' Parse coordinate filename to extract metadata
#' @noRd
parse_coord_filename <- function(path) {
  fname <- basename(path)
  m <- regmatches(fname, regexec("^(.+)_global_coords_offset([+-]?[0-9]+)_type([12])\\.tsv\\.gz$", fname))[[1]]
  if (length(m) != 4) stop("Invalid coordinate filename: ", fname)
  list(
    organism = m[2],
    offset = as.integer(m[3]),
    type = paste0("type", m[4])
  )
}

#' Load tRNA group lookup table
#'
#' Internal function to load the pre-computed mapping of tRNA IDs to their
#' coordinate groups (offset and type).
#'
#' @param organism Character string specifying the organism.
#' @return Tibble with columns: `trna_id`, `offset`, `type`
#' @noRd
load_trna_groups <- function(organism) {
  groups_file <- system.file(
    "extdata", "coords",
    paste0(organism, "_trna_groups.tsv.gz"),
    package = "clover",
    mustWork = TRUE
  )

  readr::read_tsv(
    groups_file,
    col_types = readr::cols(
      trna_id = readr::col_character(),
      offset = readr::col_integer(),
      type = readr::col_character()
    )
  )
}

#' Load global tRNA coordinates
#'
#' Load pre-computed global coordinate mappings for tRNAs. These coordinates
#' enable alignment of tRNAs with different sequences onto a common axis
#' based on their structural positions (Sprinzl coordinates).
#'
#' tRNAs are grouped by structural type (Type I vs Type II) and labeling
#' offset to ensure proper Sprinzl position alignment. By default, all
#' groups are loaded and combined. Use `type` and `offset` parameters
#' to load specific groups.
#'
#' @param organism Character string specifying the organism. One of
#'   "ecoliK12" (E. coli K12), "sacCer" (S. cerevisiae), or "hg38" (H. sapiens).
#' @param type Optional. Filter to specific tRNA type: "type1" (standard) or
#'   "type2" (extended variable arm - Leu, Ser, Tyr).
#' @param offset Optional. Filter to specific labeling offset (integer).
#' @param path Optional path to a custom coordinate file. If provided,
#'   other parameters are ignored.
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
#'     \item{offset}{Labeling offset for this coordinate group}
#'     \item{type}{Structural type ("type1" or "type2")}
#'   }
#'
#' @export
#'
#' @examples
#' # Load all coordinate groups for yeast
#' coords <- load_global_coords("sacCer")
#'
#' # Load only Type I tRNAs
#' coords_type1 <- load_global_coords("sacCer", type = "type1")
#'
#' # Load specific group
#' coords_grp <- load_global_coords("sacCer", type = "type1", offset = 0)
load_global_coords <- function(organism = "sacCer", type = NULL, offset = NULL,
                               path = NULL) {
  if (!is.null(path)) {
    return(read_coords_file(path))
  }

  organism <- match.arg(organism, available_organisms())

  # Find matching coordinate files
  coord_dir <- system.file("extdata", "coords", package = "clover", mustWork = TRUE)

  if (!is.null(type) && !is.null(offset)) {
    # Load specific file
    offset_str <- ifelse(offset >= 0, paste0("+", offset), as.character(offset))
    # Try both formats for offset 0
    if (offset == 0) offset_str <- "0"
    pattern <- paste0("^", organism, "_global_coords_offset", offset_str, "_", type, "\\.tsv\\.gz$")
    files <- list.files(coord_dir, pattern = pattern, full.names = TRUE)
    if (length(files) == 0) {
      # Try alternate format
      offset_str <- as.character(offset)
      pattern <- paste0("^", organism, "_global_coords_offset", offset_str, "_", type, "\\.tsv\\.gz$")
      files <- list.files(coord_dir, pattern = pattern, full.names = TRUE)
    }
  } else if (!is.null(type)) {
    # Load all files for this type
    pattern <- paste0("^", organism, "_global_coords_offset[+-]?[0-9]+_", type, "\\.tsv\\.gz$")
    files <- list.files(coord_dir, pattern = pattern, full.names = TRUE)
  } else if (!is.null(offset)) {
    # Load all files for this offset
    offset_str <- ifelse(offset >= 0, paste0("[+]?", offset), as.character(offset))
    if (offset == 0) offset_str <- "[+-]?0"
    pattern <- paste0("^", organism, "_global_coords_offset", offset_str, "_type[12]\\.tsv\\.gz$")
    files <- list.files(coord_dir, pattern = pattern, full.names = TRUE)
  } else {
    # Load all files for organism
    pattern <- paste0("^", organism, "_global_coords_offset[+-]?[0-9]+_type[12]\\.tsv\\.gz$")
    files <- list.files(coord_dir, pattern = pattern, full.names = TRUE)
  }

  if (length(files) == 0) {
    stop("No coordinate files found matching criteria for organism: ", organism)
  }

  # Load and combine all matching files
  purrr::map_dfr(files, function(f) {
    meta <- parse_coord_filename(f)
    read_coords_file(f) |>
      dplyr::mutate(
        offset = meta$offset,
        type = meta$type
      )
  })
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
#' Automatically classifies each tRNA into the correct coordinate group
#' based on pre-computed lookup tables.
#'
#' @param bcerror Tibble of bcerror data from [read_bcerror()]
#' @param organism Character string specifying the organism name.
#'
#' @return The input bcerror tibble with additional columns:
#'   `global_index`, `sprinzl_label`, `region`, `offset`, `type`, and `is_adapter`.
#'
#' @details
#' **Auto-Classification**:
#'
#' Each tRNA is automatically matched to the appropriate coordinate group

#' (offset + type combination) using pre-computed lookup tables. tRNAs not
#' found in any group (e.g., SeC, mitochondrial, initiator Met) will have
#' NA values for coordinate columns.
#'
#' **Adapter Offset Handling**:
#'
#' The join accounts for a 24-nucleotide adapter sequence present at the
#' 5' end of reference sequences used in the nanopore tRNA-seq protocol.
#'
#' Reference FASTA structure: `[24nt adapter][73nt mature tRNA][33nt tail]`
#'
#' - **bcerror positions**: 1-based from start of full reference (adapter + tRNA + tail)
#' - **coordinate seq_index**: 1-based from start of mature tRNA sequence
#' - **Offset applied**: seq_index + 24 matches bcerror position
#'
#' @export
#'
#' @examples
#' bcerr_path <- clover_example("yeast/grande.bcerr.tsv.gz")
#' bcerr <- read_bcerror(bcerr_path)
#' bcerr_with_coords <- add_global_coords(bcerr, "sacCer")
#' head(bcerr_with_coords)
add_global_coords <- function(bcerror, organism) {
  organism <- match.arg(organism, available_organisms())

  # Load lookup table for tRNA -> group mapping
  trna_groups <- load_trna_groups(organism)

  # Find which groups are needed for this data
  bcerror_trnas <- unique(bcerror$ref)
  matched_groups <- trna_groups |>
    dplyr::filter(trna_id %in% bcerror_trnas) |>
    dplyr::distinct(offset, type)

  # Warn about unmatched tRNAs
  unmatched <- setdiff(bcerror_trnas, trna_groups$trna_id)
  if (length(unmatched) > 0) {
    warning(
      length(unmatched), " tRNA(s) not found in coordinate groups ",
      "(may be SeC, mitochondrial, or initiator Met): ",
      paste(head(unmatched, 3), collapse = ", "),
      if (length(unmatched) > 3) paste0(", ... (", length(unmatched) - 3, " more)")
    )
  }

  if (nrow(matched_groups) == 0) {
    warning("No tRNAs matched any coordinate groups")
    return(
      bcerror |>
        dplyr::mutate(
          global_index = NA_integer_,
          sprinzl_label = NA_character_,
          region = NA_character_,
          offset = NA_integer_,
          type = NA_character_,
          is_adapter = TRUE
        )
    )
  }

  # Load only the needed coordinate files
  coord_dir <- system.file("extdata", "coords", package = "clover", mustWork = TRUE)
  coords_list <- purrr::map2(matched_groups$offset, matched_groups$type, function(off, tp) {
    offset_str <- as.character(off)
    if (off > 0) offset_str <- paste0("+", off)
    pattern <- paste0("^", organism, "_global_coords_offset", offset_str, "_", tp, "\\.tsv\\.gz$")
    files <- list.files(coord_dir, pattern = pattern, full.names = TRUE)
    if (length(files) == 0) {
      # Try without + for positive offsets
      offset_str <- as.character(off)
      pattern <- paste0("^", organism, "_global_coords_offset", offset_str, "_", tp, "\\.tsv\\.gz$")
      files <- list.files(coord_dir, pattern = pattern, full.names = TRUE)
    }
    if (length(files) == 0) return(NULL)
    read_coords_file(files[1]) |>
      dplyr::mutate(offset = off, type = tp)
  })

  coords <- dplyr::bind_rows(coords_list)

  # Adjust for 24nt adapter: bcerror pos includes adapter, coords seq_index does not
  coords_adj <- coords |>
    dplyr::mutate(seq_index_with_adapter = seq_index + 24L)

  result <- dplyr::left_join(
    bcerror,
    coords_adj,
    by = c("ref" = "trna_id", "pos" = "seq_index_with_adapter")
  ) |>
    dplyr::select(-seq_index)  # Remove seq_index, keep only pos

  # Mark adapter regions explicitly for filtering in downstream analysis
  # 5' adapter: positions 1-24
  # 3' tail: positions beyond mature tRNA (where global_index is NA)
  result <- result |>
    dplyr::mutate(
      is_adapter = pos <= 24 | is.na(global_index)
    )

  result
}

#' Get unique global index labels for plotting
#'
#' Returns a mapping of global indices to Sprinzl labels for axis annotation.
#' When coordinates contain multiple groups, returns labels from each group
#' separately since global indices have different meanings per group.
#'
#' @param coords Tibble of global coordinates from [load_global_coords()]
#'
#' @return A named character vector where names are global indices and
#'   values are Sprinzl labels. If multiple groups present, returns a list
#'   of named vectors keyed by "offset_type".
#'
#' @export
#'
#' @examples
#' coords <- load_global_coords("sacCer", type = "type1", offset = 0)
#' labels <- get_global_labels(coords)
#' head(labels)
get_global_labels <- function(coords) {
  # Check if coords has multiple groups
  if ("offset" %in% names(coords) && "type" %in% names(coords)) {
    groups <- coords |>
      dplyr::distinct(offset, type)

    if (nrow(groups) > 1) {
      # Return list of label vectors by group
      labels_list <- lapply(seq_len(nrow(groups)), function(i) {
        grp_coords <- coords |>
          dplyr::filter(offset == groups$offset[i], type == groups$type[i])
        label_map <- grp_coords |>
          dplyr::filter(!is.na(sprinzl_label), sprinzl_label != "-1", sprinzl_label != "") |>
          dplyr::distinct(global_index, sprinzl_label) |>
          dplyr::arrange(global_index)
        stats::setNames(label_map$sprinzl_label, label_map$global_index)
      })
      names(labels_list) <- paste0("offset", groups$offset, "_", groups$type)
      return(labels_list)
    }
  }

  # Single group - return simple named vector
  label_map <- coords |>
    dplyr::filter(!is.na(sprinzl_label), sprinzl_label != "-1", sprinzl_label != "") |>
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
#' coords <- load_global_coords("sacCer", type = "type1", offset = 0)
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


#' Classify tRNA as Type I or Type II
#'
#' Type II tRNAs have extended variable loops (9-24 extra nucleotides). These
#' include Leucine, Serine, and Tyrosine tRNAs. All other nuclear tRNAs are
#' Type I with standard 73-77 nucleotide length.
#'
#' Note: Selenocysteine (SeC), mitochondrial, and initiator Met tRNAs are
#' excluded from the coordinate system entirely due to structural incompatibility.
#'
#' @param trna_ids Character vector of tRNA IDs (e.g., "tRNA-Leu-CAA-1-1")
#' @return Character vector of "Type I" or "Type II"
#'
#' @export
#'
#' @examples
#' # Type II tRNAs
#' classify_trna_type(c("tRNA-Leu-CAA-1-1", "tRNA-Ser-GCT-1-1"))
#'
#' # Type I tRNAs
#' classify_trna_type(c("tRNA-Ala-GGC-1-1", "tRNA-Arg-ACG-1-1"))
classify_trna_type <- function(trna_ids) {
  ifelse(
    grepl("tRNA-(Leu|Ser|Tyr)", trna_ids),
    "Type II",
    "Type I"
  )
}
