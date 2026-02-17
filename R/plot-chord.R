# Circlize chord diagram functions ---------------------------------------------

#' Plot a chord diagram of modification co-occurrence.
#'
#' Display a circlize chord diagram showing significant pairwise modification
#' co-occurrence (odds ratios) for a single sample. Chords are colored by
#' direction: positive odds ratios (co-occurring modifications) vs. negative
#' (mutually exclusive). Chord width reflects the magnitude of the log odds
#' ratio.
#'
#' @param odds_data A tibble of odds ratio data with columns `pos1`, `pos2`,
#'   `log_odds_ratio` (or column specified by `or_col`), `p_value` (or column
#'   specified by `p_col`), and `total_obs`.
#' @param or_col Column name (string) for the odds ratio value.
#'   Default `"log_odds_ratio"`.
#' @param or_cutoff Minimum absolute value of `or_col` for a chord to be
#'   drawn. Default `1.0`.
#' @param p_col Column name (string) for the p-value. Default `"p_value"`.
#' @param p_cutoff Maximum p-value for a chord to be drawn. Default `0.01`.
#' @param min_obs Minimum number of observations (`total_obs`) for a pair
#'   to be included. Default `100`.
#' @param positive_color Color for positive odds ratios (co-occurrence).
#'   Default `"#D55E00"` (vermillion).
#' @param negative_color Color for negative odds ratios (exclusion).
#'   Default `"#0072B2"` (blue).
#' @param sprinzl_coords An optional tibble from [read_sprinzl_coords()] used
#'   to order sectors by Sprinzl position and color by structural region.
#' @param mods An optional tibble of modification annotations (e.g., from
#'   [fetch_modomics_mods()]) with columns `pos` (seq_index) and `mod1`.
#'   When provided along with `sprinzl_coords`, modification positions are
#'   highlighted as an annotation ring.
#' @param title Optional plot title.
#' @param transparency Transparency for chord colors (0 = opaque, 1 = fully
#'   transparent). Default `0.4`.
#'
#' @return `invisible(NULL)`. Called for its side effect of producing a
#'   base graphics plot.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' or_data <- read_odds_ratios("sample1.odds_ratios.tsv.gz")
#' plot_chord_or(or_data)
#' }
plot_chord_or <- function(
  odds_data,
  or_col = "log_odds_ratio",
  or_cutoff = 1.0,
  p_col = "p_value",
  p_cutoff = 0.01,
  min_obs = 100,
  positive_color = "#D55E00",
  negative_color = "#0072B2",
  sprinzl_coords = NULL,
  mods = NULL,
  title = NULL,
  transparency = 0.4
) {
  rlang::check_installed("circlize", reason = "to create chord diagrams.")

  # Filter significant pairs (drop Inf values that break circlize)
  sig_data <- odds_data |>
    dplyr::filter(
      is.finite(.data[[or_col]]),
      abs(.data[[or_col]]) >= or_cutoff,
      .data[[p_col]] <= p_cutoff,
      total_obs >= min_obs
    )

  if (nrow(sig_data) == 0) {
    cli_inform("No significant pairs found with current cutoffs.")
    return(invisible(NULL))
  }

  # Prepare chord data frame
  chord_df <- data.frame(
    from = as.character(sig_data$pos1),
    to = as.character(sig_data$pos2),
    value = abs(sig_data[[or_col]]),
    stringsAsFactors = FALSE
  )

  # Set chord colors based on sign
  chord_colors <- ifelse(
    sig_data[[or_col]] > 0,
    positive_color,
    negative_color
  )

  # Map to Sprinzl labels if coords provided
  if (!is.null(sprinzl_coords)) {
    chord_df <- .map_to_sprinzl(chord_df, sprinzl_coords)
    # Update chord_colors to match filtered rows
    chord_colors <- chord_colors[seq_len(nrow(chord_df))]
  }

  # Set up sectors and grid colors
  setup <- .setup_chord_sectors(chord_df, sprinzl_coords)

  # Build adjacency matrix so all sectors appear (even without chords)
  adj_mat <- .build_adjacency_matrix(chord_df, setup$order)

  # Equalize sector widths with diagonal padding
  eq <- .equalize_sectors(adj_mat)

  # Map chord colors to match adjacency matrix links
  col_mat <- .build_color_matrix(
    chord_df,
    chord_colors,
    setup$order,
    transparency
  )
  # Pad color matrix diagonal to match equalized adjacency matrix
  diag(col_mat) <- grDevices::rgb(1, 1, 1, alpha = 0)

  # Draw chord diagram
  circlize::circos.clear()
  circlize::circos.par(start.degree = 90, gap.degree = 2)

  circlize::chordDiagram(
    eq$mat,
    order = setup$order,
    grid.col = setup$grid_col,
    col = col_mat,
    transparency = 0,
    annotationTrack = "grid",
    preAllocateTracks = list(track.height = 0.05),
    reduce = -1,
    self.link = 1,
    link.visible = eq$link_visible
  )

  # Add sector labels
  .add_sector_labels()

  # Add chord color legend
  .add_chord_legend(
    labels = c("Co-occurring", "Exclusive"),
    colors = c(positive_color, negative_color)
  )

  # Add annotation rings when sprinzl coords provided
  if (!is.null(sprinzl_coords)) {
    .add_region_ring(setup$regions)
    .add_nucleotide_ring(setup$residues)

    if (!is.null(mods)) {
      .add_modification_ring(mods, sprinzl_coords, setup$order)
    }
  }

  if (!is.null(title)) {
    graphics::title(title)
  }

  circlize::circos.clear()
  invisible(NULL)
}

#' Compute ratio of odds ratios between conditions.
#'
#' Compare modification co-occurrence between two conditions by computing
#' the ratio of odds ratios (ROR). Replicates within each condition are
#' aggregated using the specified function.
#'
#' @param odds_data A combined tibble of odds ratio data with a `condition`
#'   column (or column specified by `condition_col`) and `sample_id`.
#' @param condition_col Column name (string) for condition labels.
#'   Default `"condition"`.
#' @param numerator Value of `condition_col` for the numerator condition.
#' @param denominator Value of `condition_col` for the denominator condition.
#' @param min_obs Minimum `total_obs` for a pair to be included.
#'   Default `100`.
#' @param agg_fun Function to aggregate replicate log odds ratios.
#'   Default `mean`.
#'
#' @return A tibble with columns: `pos1`, `pos2`, `or_numerator`,
#'   `or_denominator`, `ror`, and `log_ror`.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' combined_or <- read_odds_ratios_multi(paths)
#' combined_or$condition <- ifelse(
#'   grepl("wt", combined_or$sample_id), "wt", "mut"
#' )
#' ror <- compute_ror(combined_or, numerator = "mut", denominator = "wt")
#' }
compute_ror <- function(
  odds_data,
  condition_col = "condition",
  numerator,
  denominator,
  min_obs = 100,
  agg_fun = mean
) {
  # Filter by minimum observations
  filtered <- odds_data |>
    dplyr::filter(total_obs >= min_obs)

  # Aggregate replicates within each condition
  agg <- filtered |>
    dplyr::group_by(
      .data[[condition_col]],
      pos1,
      pos2
    ) |>
    dplyr::summarise(
      mean_log_or = agg_fun(log_odds_ratio),
      .groups = "drop"
    )

  # Separate numerator and denominator
  num <- agg |>
    dplyr::filter(.data[[condition_col]] == numerator) |>
    dplyr::select(pos1, pos2, or_numerator = mean_log_or)

  denom <- agg |>
    dplyr::filter(.data[[condition_col]] == denominator) |>
    dplyr::select(pos1, pos2, or_denominator = mean_log_or)

  # Join and compute ROR
  dplyr::inner_join(num, denom, by = c("pos1", "pos2")) |>
    dplyr::mutate(
      log_ror = or_numerator - or_denominator,
      ror = exp(log_ror)
    )
}

#' Plot a chord diagram of modification rewiring between conditions.
#'
#' Display a circlize chord diagram showing position pairs where modification
#' co-occurrence has changed between two conditions, based on the ratio of
#' odds ratios (ROR). Chords are colored by direction: gained dependencies
#' (positive log ROR) vs. lost dependencies (negative log ROR).
#'
#' @param ror_data A tibble from [compute_ror()] with columns `pos1`, `pos2`,
#'   and `log_ror`.
#' @param ror_cutoff Minimum absolute value of `log_ror` for a chord to be
#'   drawn. Default `0.5`.
#' @param gained_color Color for gained dependencies (positive log ROR).
#'   Default `"#D55E00"` (vermillion).
#' @param lost_color Color for lost dependencies (negative log ROR).
#'   Default `"#0072B2"` (blue).
#' @param sprinzl_coords An optional tibble from [read_sprinzl_coords()] used
#'   to order sectors and color by structural region.
#' @param mods An optional tibble of modification annotations (e.g., from
#'   [fetch_modomics_mods()]) with columns `pos` (seq_index) and `mod1`.
#'   When provided along with `sprinzl_coords`, modification positions are
#'   highlighted as an annotation ring.
#' @param title Optional plot title.
#' @param transparency Transparency for chord colors. Default `0.4`.
#'
#' @return `invisible(NULL)`. Called for its side effect of producing a
#'   base graphics plot.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' ror <- compute_ror(combined_or, numerator = "mut", denominator = "wt")
#' plot_chord_ror(ror)
#' }
plot_chord_ror <- function(
  ror_data,
  ror_cutoff = 0.5,
  gained_color = "#D55E00",
  lost_color = "#0072B2",
  sprinzl_coords = NULL,
  mods = NULL,
  title = NULL,
  transparency = 0.4
) {
  rlang::check_installed("circlize", reason = "to create chord diagrams.")

  # Filter by ROR cutoff (drop Inf values that break circlize)
  sig_data <- ror_data |>
    dplyr::filter(is.finite(log_ror), abs(log_ror) >= ror_cutoff)

  if (nrow(sig_data) == 0) {
    cli_inform("No pairs exceed the ROR cutoff.")
    return(invisible(NULL))
  }

  # Prepare chord data frame
  chord_df <- data.frame(
    from = as.character(sig_data$pos1),
    to = as.character(sig_data$pos2),
    value = abs(sig_data$log_ror),
    stringsAsFactors = FALSE
  )

  # Set chord colors based on direction
  chord_colors <- ifelse(
    sig_data$log_ror > 0,
    gained_color,
    lost_color
  )

  # Map to Sprinzl labels if coords provided
  if (!is.null(sprinzl_coords)) {
    chord_df <- .map_to_sprinzl(chord_df, sprinzl_coords)
    chord_colors <- chord_colors[seq_len(nrow(chord_df))]
  }

  # Set up sectors and grid colors
  setup <- .setup_chord_sectors(chord_df, sprinzl_coords)

  # Build adjacency matrix so all sectors appear (even without chords)
  adj_mat <- .build_adjacency_matrix(chord_df, setup$order)

  # Equalize sector widths with diagonal padding
  eq <- .equalize_sectors(adj_mat)

  # Map chord colors to match adjacency matrix links
  col_mat <- .build_color_matrix(
    chord_df,
    chord_colors,
    setup$order,
    transparency
  )
  # Pad color matrix diagonal to match equalized adjacency matrix
  diag(col_mat) <- grDevices::rgb(1, 1, 1, alpha = 0)

  # Draw chord diagram
  circlize::circos.clear()
  circlize::circos.par(start.degree = 90, gap.degree = 2)

  circlize::chordDiagram(
    eq$mat,
    order = setup$order,
    grid.col = setup$grid_col,
    col = col_mat,
    transparency = 0,
    annotationTrack = "grid",
    preAllocateTracks = list(track.height = 0.05),
    reduce = -1,
    self.link = 1,
    link.visible = eq$link_visible
  )

  # Add sector labels
  .add_sector_labels()

  # Add chord color legend
  .add_chord_legend(
    labels = c("Gained", "Lost"),
    colors = c(gained_color, lost_color)
  )

  # Add annotation rings when sprinzl coords provided
  if (!is.null(sprinzl_coords)) {
    .add_region_ring(setup$regions)
    .add_nucleotide_ring(setup$residues)

    if (!is.null(mods)) {
      .add_modification_ring(mods, sprinzl_coords, setup$order)
    }
  }

  if (!is.null(title)) {
    graphics::title(title)
  }

  circlize::circos.clear()
  invisible(NULL)
}

# Internal helpers -----------------------------------------------------------

#' Map pos1/pos2 from seq_index to Sprinzl labels.
#'
#' @param chord_df Data frame with `from` and `to` columns (seq_index as
#'   character).
#' @param sprinzl_coords Tibble from [read_sprinzl_coords()].
#'
#' @return `chord_df` with `from`/`to` replaced by Sprinzl labels. Rows
#'   where either position has no Sprinzl mapping are dropped.
#' @noRd
.map_to_sprinzl <- function(chord_df, sprinzl_coords) {
  lookup <- sprinzl_coords |>
    dplyr::select(seq_index, sprinzl_label) |>
    dplyr::distinct(seq_index, .keep_all = TRUE) |>
    dplyr::filter(!is.na(sprinzl_label))

  chord_df$from_idx <- as.numeric(chord_df$from)
  chord_df$to_idx <- as.numeric(chord_df$to)

  chord_df <- chord_df |>
    dplyr::left_join(lookup, by = c("from_idx" = "seq_index")) |>
    dplyr::rename(from_label = sprinzl_label) |>
    dplyr::left_join(lookup, by = c("to_idx" = "seq_index")) |>
    dplyr::rename(to_label = sprinzl_label)

  # Drop pairs where either position has no mapping

  chord_df <- chord_df |>
    dplyr::filter(!is.na(from_label), !is.na(to_label))

  chord_df$from <- chord_df$from_label
  chord_df$to <- chord_df$to_label

  chord_df |>
    dplyr::select(-from_idx, -to_idx, -from_label, -to_label)
}

#' Set up chord diagram sectors from position pairs.
#' @noRd
.setup_chord_sectors <- function(chord_df, sprinzl_coords = NULL) {
  if (!is.null(sprinzl_coords)) {
    # Use ALL non-NA Sprinzl positions as sectors
    all_labels <- sprinzl_coords |>
      dplyr::filter(!is.na(sprinzl_label)) |>
      dplyr::pull(sprinzl_label) |>
      unique()

    ordered <- order_sprinzl_positions(all_labels)
    sector_order <- levels(ordered)

    # Build region mapping for the outer ring
    region_map <- sprinzl_coords |>
      dplyr::filter(!is.na(sprinzl_label)) |>
      dplyr::select(sprinzl_label, region) |>
      dplyr::distinct(sprinzl_label, .keep_all = TRUE)

    regions <- stats::setNames(region_map$region, region_map$sprinzl_label)

    # Build residue mapping for the nucleotide ring
    residue_map <- sprinzl_coords |>
      dplyr::filter(!is.na(sprinzl_label)) |>
      dplyr::select(sprinzl_label, residue) |>
      dplyr::distinct(sprinzl_label, .keep_all = TRUE)

    residues <- stats::setNames(
      residue_map$residue,
      residue_map$sprinzl_label
    )

    # Neutral grid colors; region coloring goes on outer ring
    grid_col <- rep("grey90", length(sector_order))
    names(grid_col) <- sector_order
  } else {
    all_positions <- unique(c(chord_df$from, chord_df$to))
    # Sort positions numerically where possible
    sector_order <- all_positions[order(
      suppressWarnings(as.numeric(all_positions)),
      all_positions
    )]
    grid_col <- rep("grey70", length(sector_order))
    names(grid_col) <- sector_order
    regions <- NULL
    residues <- NULL
  }

  list(
    order = sector_order,
    grid_col = grid_col,
    regions = regions,
    residues = residues
  )
}

#' Build adjacency matrix from chord data.
#'
#' Creates a square matrix with all sector positions as rows/columns.
#' Positions without connections have all-zero rows/columns but still
#' appear as sectors in the diagram.
#' @noRd
.build_adjacency_matrix <- function(chord_df, sector_order) {
  n <- length(sector_order)
  mat <- matrix(0, nrow = n, ncol = n)
  rownames(mat) <- sector_order
  colnames(mat) <- sector_order

  for (i in seq_len(nrow(chord_df))) {
    ri <- match(chord_df$from[i], sector_order)
    ci <- match(chord_df$to[i], sector_order)
    if (!is.na(ri) && !is.na(ci)) {
      mat[ri, ci] <- chord_df$value[i]
    }
  }

  mat
}

#' Build color matrix for adjacency matrix chord diagram.
#'
#' Maps per-chord colors into a matrix matching the adjacency matrix,
#' with transparency already applied.
#' @noRd
.build_color_matrix <- function(
  chord_df,
  chord_colors,
  sector_order,
  transparency
) {
  n <- length(sector_order)
  # Default: fully transparent (no chord)
  col_mat <- matrix(
    grDevices::rgb(1, 1, 1, alpha = 0),
    nrow = n,
    ncol = n
  )
  rownames(col_mat) <- sector_order
  colnames(col_mat) <- sector_order

  for (i in seq_len(nrow(chord_df))) {
    ri <- match(chord_df$from[i], sector_order)
    ci <- match(chord_df$to[i], sector_order)
    if (!is.na(ri) && !is.na(ci)) {
      base_col <- grDevices::col2rgb(chord_colors[i]) / 255
      col_mat[ri, ci] <- grDevices::rgb(
        base_col[1],
        base_col[2],
        base_col[3],
        alpha = 1 - transparency
      )
    }
  }

  col_mat
}

#' Add sector labels to chord diagram.
#' @noRd
.add_sector_labels <- function() {
  circlize::circos.track(
    track.index = 1,
    ylim = c(0, 1),
    panel.fun = function(x, y) {
      sector_name <- circlize::get.cell.meta.data("sector.index")
      xlim <- circlize::get.cell.meta.data("xlim")
      ylim <- circlize::get.cell.meta.data("ylim")
      circlize::circos.text(
        mean(xlim),
        ylim[1] + 0.1,
        sector_name,
        facing = "clockwise",
        niceFacing = TRUE,
        adj = c(0, 0.5),
        cex = 0.6
      )
    },
    bg.border = NA
  )
}

#' Add outer ring colored by structural region.
#' @noRd
.add_region_ring <- function(regions) {
  region_palette <- .region_colors()

  circlize::circos.track(
    ylim = c(0, 1),
    track.height = 0.05,
    bg.border = NA,
    panel.fun = function(x, y) {
      sector_name <- circlize::get.cell.meta.data("sector.index")
      xlim <- circlize::get.cell.meta.data("xlim")
      ylim <- circlize::get.cell.meta.data("ylim")

      rgn <- regions[sector_name]
      col <- if (is.na(rgn)) "grey70" else region_palette[rgn]

      circlize::circos.rect(
        xlim[1],
        ylim[1],
        xlim[2],
        ylim[2],
        col = col,
        border = NA
      )
    }
  )

  # Add legend for regions present in the data
  present_regions <- unique(stats::na.omit(regions))
  legend_colors <- region_palette[present_regions]

  # Deduplicate colors (e.g., acceptor-stem and acceptor-tail share a color)
  unique_colors <- !duplicated(legend_colors)
  legend_labels <- present_regions[unique_colors]
  legend_cols <- legend_colors[unique_colors]

  graphics::legend(
    "bottomleft",
    legend = legend_labels,
    fill = legend_cols,
    border = NA,
    bty = "n",
    cex = 0.7
  )
}

#' Add ring colored by reference nucleotide.
#' @noRd
.add_nucleotide_ring <- function(residues) {
  nuc_palette <- .nucleotide_colors()

  circlize::circos.track(
    ylim = c(0, 1),
    track.height = 0.05,
    bg.border = NA,
    panel.fun = function(x, y) {
      sector_name <- circlize::get.cell.meta.data("sector.index")
      xlim <- circlize::get.cell.meta.data("xlim")
      ylim <- circlize::get.cell.meta.data("ylim")

      nuc <- residues[sector_name]
      col <- if (is.na(nuc)) "grey90" else nuc_palette[nuc]

      circlize::circos.rect(
        xlim[1],
        ylim[1],
        xlim[2],
        ylim[2],
        col = col,
        border = NA
      )
    }
  )
}

#' Add ring highlighting modification positions.
#' @noRd
.add_modification_ring <- function(mods, sprinzl_coords, sector_order) {
  # Map mod positions (seq_index) to Sprinzl labels
  lookup <- sprinzl_coords |>
    dplyr::select(seq_index, sprinzl_label) |>
    dplyr::distinct(seq_index, .keep_all = TRUE) |>
    dplyr::filter(!is.na(sprinzl_label))

  mod_mapped <- mods |>
    dplyr::inner_join(lookup, by = c("pos" = "seq_index"))

  mod_labels <- unique(mod_mapped$sprinzl_label)

  circlize::circos.track(
    ylim = c(0, 1),
    track.height = 0.05,
    bg.border = NA,
    panel.fun = function(x, y) {
      sector_name <- circlize::get.cell.meta.data("sector.index")
      xlim <- circlize::get.cell.meta.data("xlim")
      ylim <- circlize::get.cell.meta.data("ylim")

      if (sector_name %in% mod_labels) {
        circlize::circos.rect(
          xlim[1],
          ylim[1],
          xlim[2],
          ylim[2],
          col = "#E41A1C",
          border = NA
        )
      }
    }
  )
}

#' Equalize sector widths by adding diagonal padding.
#'
#' Computes per-sector diagonal values so that every sector occupies the
#' same total width in the chord diagram (off-diagonal links + diagonal
#' padding = constant).
#'
#' @param mat Square adjacency matrix from `.build_adjacency_matrix()`.
#' @return A list with `mat` (padded matrix) and `link_visible` (logical
#'   matrix; `FALSE` on diagonal entries used for padding).
#' @noRd
.equalize_sectors <- function(mat) {
  n <- nrow(mat)
  link_total <- rowSums(mat) + colSums(mat)
  # Diagonal entries contribute to both row and column sums, so each
  # adds 2 * diag[i] to sector width. Solve for equal total widths:
  # link_total[i] + 2 * diag[i] = target
  target <- max(link_total) + 2
  diag_pad <- (target - link_total) / 2

  diag(mat) <- diag_pad

  link_visible <- matrix(TRUE, nrow = n, ncol = n)
  rownames(link_visible) <- rownames(mat)
  colnames(link_visible) <- colnames(mat)
  diag(link_visible) <- FALSE

  list(mat = mat, link_visible = link_visible)
}

#' Add a legend for chord colors.
#' @param labels Character vector of legend labels.
#' @param colors Character vector of colors matching `labels`.
#' @noRd
.add_chord_legend <- function(labels, colors) {
  graphics::legend(
    "bottomright",
    legend = labels,
    fill = colors,
    border = NA,
    bty = "n",
    cex = 0.7
  )
}

#' Named color palette for nucleotides.
#' @noRd
.nucleotide_colors <- function() {
  c(
    "A" = "#4DAF4A",
    "C" = "#377EB8",
    "G" = "#FFD92F",
    "U" = "#E41A1C"
  )
}

#' Named color palette for tRNA structural regions.
#' @noRd
.region_colors <- function() {
  c(
    "acceptor-stem" = "#E41A1C",
    "acceptor-tail" = "#E41A1C",
    "D-stem" = "#377EB8",
    "D-loop" = "#4DAF4A",
    "anticodon-stem" = "#984EA3",
    "anticodon-loop" = "#FF7F00",
    "variable-region" = "#A65628",
    "variable-arm" = "#A65628",
    "T-stem" = "#F781BF",
    "T-loop" = "#999999",
    "unknown" = "grey70"
  )
}
