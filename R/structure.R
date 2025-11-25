# Structure visualization functions -----------------------------------------------

#' Parse R2DT JSON file to extract structure coordinates
#'
#' @param json_path Path to R2DT enriched JSON file
#' @return List with `residues` (tibble with x, y, sprinzl_label) and
#'   `basepairs` (tibble with pos1, pos2, type)
#' @noRd
parse_r2dt_json <- function(json_path) {
  json <- jsonlite::read_json(json_path)


  # Extract sequence/residue data

  sequence <- json$rnaComplexes[[1]]$rnaMolecules[[1]]$sequence
  residues <- purrr::map_dfr(sequence, function(res) {
    tibble::tibble(
      residue_index = res$residueIndex,
      residue_name = res$residueName,
      x = res$x,
      y = res$y,
      sprinzl_index = res$info$templateResidueIndex %||% NA_integer_,
      sprinzl_label = res$info$templateNumberingLabel %||% NA_character_
    )
  })

  # Filter out 5' and 3' markers (residue_index 0 and max+1)
  residues <- residues |>
    dplyr::filter(
      !residue_name %in% c("5'", "3'"),
      residue_index > 0
    )

  # Extract base pairs

  basepairs_raw <- json$rnaComplexes[[1]]$rnaMolecules[[1]]$basePairs
  basepairs <- purrr::map_dfr(basepairs_raw, function(bp) {
    tibble::tibble(
      pos1 = bp$residueIndex1,
      pos2 = bp$residueIndex2,
      type = bp$basePairType %||% "canonical"
    )
  })

  list(residues = residues, basepairs = basepairs)
}


#' Load tRNA consensus structure template
#'
#' Loads pre-computed x,y coordinates for consensus tRNA secondary structure.
#' Uses a Type I tRNA template (standard 76nt) suitable for most tRNAs.
#'
#' @return List with `residues` (tibble with sprinzl_label, x, y) and
#'   `basepairs` (tibble with pos1, pos2)
#'
#' @export
#'
#' @examples
#' template <- load_structure_template()
#' head(template$residues)
load_structure_template <- function() {
  rds_path <- system.file(
    "extdata", "structure", "consensus_template.rds",
    package = "clover",
    mustWork = FALSE
  )


  if (file.exists(rds_path)) {
    return(readRDS(rds_path))
  }


  # Fall back to parsing JSON directly from bundled file
  json_path <- system.file(
    "extdata", "structure", "tRNA-Ala-GGC-1-1-B_Ala.enriched.json",
    package = "clover",
    mustWork = FALSE
  )

  if (!file.exists(json_path)) {
    stop(
      "Structure template not found. ",
      "Run create_consensus_template() to generate it."
    )
  }

  parse_r2dt_json(json_path)
}


#' Load Modomics modification reference data
#'
#' Loads known tRNA modifications from Modomics database for a given organism.
#'
#' @param organism Character string: "sacCer" (yeast), "ecoliK12", or "hg38".
#'
#' @return Tibble with columns:
#'   - `trna_id`: gtRNAdb tRNA identifier
#'   - `sprinzl_label`: Sprinzl position
#'   - `modification_short_name`: Short modification name (e.g., "D", "m1A")
#'   - `modification_name`: Full modification name
#'   - `region`: Structural region
#'
#' @export
#'
#' @examples
#' mods <- load_modomics("sacCer")
#' head(mods)
load_modomics <- function(organism = c("sacCer", "ecoliK12", "hg38")) {
  organism <- match.arg(organism)

  file_map <- c(
    sacCer = "sacCer_modomics_to_sprinzl.tsv.gz",
    ecoliK12 = "ecoli_modomics_to_sprinzl.tsv.gz",
    hg38 = "hg38_modomics_to_sprinzl.tsv.gz"
  )

  file_path <- system.file(
    "extdata", "modomics", file_map[organism],
    package = "clover",
    mustWork = TRUE
  )

  readr::read_tsv(file_path, show_col_types = FALSE) |>
    dplyr::select(
      trna_id = gtRNAdb_trna_id,
      sprinzl_label,
      modification_short_name,
      modification_name,
      region
    )
}


#' Plot tRNA secondary structure
#'
#' Creates a cloverleaf visualization of tRNA secondary structure with
#' data mapped to color. Useful for visualizing base-calling error rates or
#' modification signals at specific positions in the tRNA structure.
#'
#' For a genome-wide view across all tRNAs, see [plot_bcerror_heatmap()].
#'
#' @param data Tibble with modification calls or bcerror data. Must have
#'   `sprinzl_label` column for mapping to structure positions. Use
#'   [add_global_coords()] to add this column to bcerror data.
#' @param trna_id Character string specifying a single tRNA to plot.
#'   If NULL, aggregates across all tRNAs in data by taking the mean
#'   value per Sprinzl position.
#' @param fill Column name to map to fill color. Default is "error_rate".
#'   Can be any numeric column in data (e.g., "mis", "ins", "del",
#'   "bcerror_residual").
#' @param compare Character vector of condition names for faceting.
#'   Data must have a `condition` column matching these values.
#' @param modomics Logical or tibble. If TRUE, loads Modomics for organism
#'   and overlays known modification positions as red circles. If tibble,
#'   uses provided Modomics data. Set to FALSE to disable overlay.
#' @param organism Organism code for Modomics lookup when `modomics = TRUE`.
#'   One of "sacCer", "ecoliK12", or "hg38".
#' @param show_labels Character: "sprinzl" to label positions, "residue" to
#'   show nucleotides, or "none" for no labels.
#' @param point_size Numeric size for residue points.
#'
#' @return A ggplot2 object displaying the tRNA cloverleaf structure. Positions
#'   with higher fill values appear brighter (viridis magma palette). Red
#'   circles indicate known modification sites when `modomics = TRUE`.
#'
#' @export
#'
#' @examples
#' # Using bundled example data
#' bcerr <- read_bcerror(clover_example("yeast/grande.bcerr.tsv.gz")) |>
#'   add_global_coords("sacCer")
#'
#' # Plot single tRNA with Modomics overlay
#' plot_trna_structure(
#'   bcerr,
#'   trna_id = "nuc-tRNA-Ala-AGC-1-1",
#'   fill = "error_rate",
#'   modomics = TRUE,
#'   organism = "sacCer"
#' )
#'
#' # Aggregate across all tRNAs (consensus view)
#' plot_trna_structure(bcerr, fill = "error_rate")
#'
#' \dontrun{
#' # Compare conditions with faceting
#' bcerr_grande <- read_bcerror("grande.bcerr.tsv.gz") |>
#'   add_global_coords("sacCer") |>
#'   dplyr::mutate(condition = "grande")
#'
#' bcerr_petite <- read_bcerror("petite.bcerr.tsv.gz") |>
#'   add_global_coords("sacCer") |>
#'   dplyr::mutate(condition = "petite")
#'
#' bcerr_both <- dplyr::bind_rows(bcerr_grande, bcerr_petite)
#'
#' plot_trna_structure(
#'   bcerr_both,
#'   compare = c("grande", "petite"),
#'   fill = "error_rate"
#' )
#' }
plot_trna_structure <- function(
    data,
    trna_id = NULL,
    fill = "error_rate",
    compare = NULL,
    modomics = TRUE,
    organism = "sacCer",
    show_labels = c("none", "sprinzl", "residue"),
    point_size = 4
) {
  show_labels <- match.arg(show_labels)

  # Validate required column

if (!"sprinzl_label" %in% names(data)) {
    stop("data must have 'sprinzl_label' column. Use add_global_coords() first.")
  }

  if (!fill %in% names(data)) {
    stop("Column '", fill, "' not found in data.")
  }

  # Load structure template
  template <- load_structure_template()

  # Filter to specific tRNA or aggregate
  if (!is.null(trna_id)) {
    ref_col <- if ("ref" %in% names(data)) "ref" else "trna_id"
    data <- dplyr::filter(data, .data[[ref_col]] == trna_id)
    if (nrow(data) == 0) {
      stop("No data found for tRNA: ", trna_id)
    }
  } else {
    # Aggregate: mean of fill value per sprinzl position
    group_cols <- "sprinzl_label"
    if (!is.null(compare) && "condition" %in% names(data)) {
      group_cols <- c(group_cols, "condition")
    }
    data <- data |>
      dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) |>
      dplyr::summarize(
        .fill_value = mean(.data[[fill]], na.rm = TRUE),
        .groups = "drop"
      )
    # Rename .fill_value back to original column name
    names(data)[names(data) == ".fill_value"] <- fill
  }

  # Join data to template coordinates
  plot_data <- template$residues |>
    dplyr::left_join(
      data,
      by = "sprinzl_label"
    )

  # Prepare base pair line coordinates
  bp_coords <- template$basepairs |>
    dplyr::left_join(
      template$residues |> dplyr::select(residue_index, x, y),
      by = c("pos1" = "residue_index")
    ) |>
    dplyr::rename(x1 = x, y1 = y) |>
    dplyr::left_join(
      template$residues |> dplyr::select(residue_index, x, y),
      by = c("pos2" = "residue_index")
    ) |>
    dplyr::rename(x2 = x, y2 = y)

  # Build plot
  p <- ggplot() +
    # Base pair lines (backbone)
    geom_segment(
      data = bp_coords,
      aes(x = x1, y = y1, xend = x2, yend = y2),
      color = "grey70",
      linewidth = 0.5
    ) +
    # Nucleotide positions
    geom_point(
      data = plot_data,
      aes(x = x, y = y, fill = .data[[fill]]),
      shape = 21,
      size = point_size,
      color = "black",
      stroke = 0.3
    ) +
    scale_fill_viridis_c(
      option = "magma",
      na.value = "grey90",
      name = fill
    ) +
    coord_fixed() +
    theme_void() +
    theme(
      legend.position = "right",
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA)
    )

  # Add Modomics overlay if requested
  if (!isFALSE(modomics)) {
    if (isTRUE(modomics)) {
      modomics_data <- load_modomics(organism)
    } else {
      modomics_data <- modomics
    }

    # Get positions with known modifications
    known_mods <- modomics_data |>
      dplyr::distinct(sprinzl_label)

    modomics_coords <- template$residues |>
      dplyr::semi_join(known_mods, by = "sprinzl_label")

    if (nrow(modomics_coords) > 0) {
      p <- p +
        geom_point(
          data = modomics_coords,
          aes(x = x, y = y),
          shape = 1,
          size = point_size + 2,
          color = "red",
          stroke = 0.8
        )
    }
  }

  # Add labels if requested
  if (show_labels == "sprinzl") {
    label_data <- plot_data |>
      dplyr::filter(as.numeric(sprinzl_label) %% 10 == 0 | sprinzl_label == "1")
    p <- p +
      geom_text(
        data = label_data,
        aes(x = x, y = y, label = sprinzl_label),
        size = 2.5,
        vjust = -1.5
      )
  } else if (show_labels == "residue") {
    p <- p +
      geom_text(
        data = plot_data,
        aes(x = x, y = y, label = residue_name),
        size = 2,
        fontface = "bold"
      )
  }

  # Add faceting for condition comparison
  if (!is.null(compare) && "condition" %in% names(data)) {
    p <- p + facet_wrap(~condition, ncol = length(compare))
  }

  p
}


#' Read modification calls from upstream pipeline
#'
#' Imports consensus modification calls from the tRNA modification detection
#' pipeline. These files contain positions called as modified in at least
#' 2 of 3 replicates.
#'
#' @param path Path to consensus TSV file (e.g., "yeast_WT_consensus.tsv")
#' @param organism Organism code for coordinate lookup: "sacCer", "ecoliK12", "hg38"
#'
#' @return Tibble with standardized columns:
#'   - `ref`: tRNA identifier
#'   - `pos`: Position in tRNA sequence
#'   - `sprinzl_label`: Sprinzl position
#'   - `region`: Structural region
#'   - `bcerror_residual`: Mean delta BCError (native - IVT)
#'   - `p_adjusted`: Mean adjusted p-value
#'   - `n_replicates`: Number of replicates calling this position
#'   - `is_consensus`: TRUE if called in >= 2 replicates
#'   - `is_modomics`: TRUE if position has known Modomics modification
#'
#' @export
#'
#' @examples
#' \dontrun{
#' mods <- read_mod_calls("yeast_WT_consensus.tsv", organism = "sacCer")
#' plot_trna_structure(mods, fill = "bcerror_residual")
#' }
read_mod_calls <- function(path, organism = c("sacCer", "ecoliK12", "hg38")) {
  organism <- match.arg(organism)

  raw <- readr::read_tsv(path, show_col_types = FALSE)

  # Standardize column names (handle variations)
  # Expected columns: Reference, Position, sprinzl_label, region,
  #                   BCErrorFreq_residual, p_adjusted, n_replicates,
  #                   is_consensus, is_modomics

  # Rename to standard names
  result <- raw |>
    dplyr::rename_with(
      ~ dplyr::case_when(
        . == "Reference" ~ "ref",
        . == "Position" ~ "pos",
        . == "BCErrorFreq_residual" ~ "bcerror_residual",
        TRUE ~ .
      )
    )

  # Add global_index from coordinates if available
  coords <- load_global_coords(organism)

  result <- result |>
    dplyr::left_join(
      coords |> dplyr::select(trna_id, seq_index, global_index),
      by = c("ref" = "trna_id", "pos" = "seq_index")
    )

  result
}


#' Create consensus template from R2DT JSON
#'
#' Helper function to generate the consensus template RDS file from an
#' R2DT JSON. Run this once to create the package data file.
#'
#' @param json_path Path to R2DT enriched JSON file
#' @param output_path Path to save RDS file
#'
#' @return Invisibly returns the template list
#' @noRd
create_consensus_template <- function(
    json_path,
    output_path = "inst/extdata/structure/consensus_template.rds"
) {
  template <- parse_r2dt_json(json_path)
  saveRDS(template, output_path)
  invisible(template)
}
