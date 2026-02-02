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
#'   drawn. Default `0.5`.
#' @param p_col Column name (string) for the p-value. Default `"p_value"`.
#' @param p_cutoff Maximum p-value for a chord to be drawn. Default `0.05`.
#' @param min_obs Minimum number of observations (`total_obs`) for a pair
#'   to be included. Default `50`.
#' @param positive_color Color for positive odds ratios (co-occurrence).
#'   Default `"#D55E00"` (vermillion).
#' @param negative_color Color for negative odds ratios (exclusion).
#'   Default `"#0072B2"` (blue).
#' @param sprinzl_coords An optional tibble from [read_sprinzl_coords()] used
#'   to order sectors by Sprinzl position and color by structural region.
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
  or_cutoff = 0.5,
  p_col = "p_value",
  p_cutoff = 0.05,
  min_obs = 50,
  positive_color = "#D55E00",
  negative_color = "#0072B2",
  sprinzl_coords = NULL,
  title = NULL,
  transparency = 0.4
) {
  rlang::check_installed("circlize", reason = "to create chord diagrams.")

  # Filter significant pairs
  sig_data <- odds_data |>
    dplyr::filter(
      abs(.data[[or_col]]) >= or_cutoff,
      .data[[p_col]] <= p_cutoff,
      total_obs >= min_obs
    )

  if (nrow(sig_data) == 0) {
    message("No significant pairs found with current cutoffs.")
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

  # Set up sectors and grid colors
  setup <- .setup_chord_sectors(chord_df, sprinzl_coords)

  # Draw chord diagram
  circlize::circos.clear()
  circlize::circos.par(start.degree = 90, gap.degree = 2)

  circlize::chordDiagram(
    chord_df,
    order = setup$order,
    grid.col = setup$grid_col,
    col = chord_colors,
    transparency = transparency,
    annotationTrack = "grid",
    preAllocateTracks = list(track.height = 0.05)
  )

  # Add sector labels
  circlize::circos.track(
    track.index = 1,
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
#'   Default `50`.
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
  min_obs = 50,
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
  title = NULL,
  transparency = 0.4
) {
  rlang::check_installed("circlize", reason = "to create chord diagrams.")

  # Filter by ROR cutoff
  sig_data <- ror_data |>
    dplyr::filter(abs(log_ror) >= ror_cutoff)

  if (nrow(sig_data) == 0) {
    message("No pairs exceed the ROR cutoff.")
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

  # Set up sectors and grid colors
  setup <- .setup_chord_sectors(chord_df, sprinzl_coords)

  # Draw chord diagram
  circlize::circos.clear()
  circlize::circos.par(start.degree = 90, gap.degree = 2)

  circlize::chordDiagram(
    chord_df,
    order = setup$order,
    grid.col = setup$grid_col,
    col = chord_colors,
    transparency = transparency,
    annotationTrack = "grid",
    preAllocateTracks = list(track.height = 0.05)
  )

  # Add sector labels
  circlize::circos.track(
    track.index = 1,
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

  if (!is.null(title)) {
    graphics::title(title)
  }

  circlize::circos.clear()
  invisible(NULL)
}

# Internal helpers -----------------------------------------------------------

#' Set up chord diagram sectors from position pairs.
#' @noRd
.setup_chord_sectors <- function(chord_df, sprinzl_coords = NULL) {
  all_positions <- unique(c(chord_df$from, chord_df$to))

  if (!is.null(sprinzl_coords)) {
    # Use Sprinzl ordering
    ordered <- order_sprinzl_positions(all_positions)
    sector_order <- levels(ordered)
    sector_order <- sector_order[sector_order %in% all_positions]

    # Color sectors by structural region
    region_map <- sprinzl_coords |>
      dplyr::select(sprinzl_label, region) |>
      dplyr::distinct(sprinzl_label, .keep_all = TRUE)

    region_palette <- .region_colors()

    grid_col <- vapply(
      sector_order,
      function(pos) {
        rgn <- region_map$region[region_map$sprinzl_label == pos]
        if (length(rgn) == 0) {
          return("grey70")
        }
        region_palette[rgn[1]]
      },
      character(1)
    )
    names(grid_col) <- sector_order
  } else {
    # Sort positions numerically where possible
    sector_order <- all_positions[order(
      suppressWarnings(as.numeric(all_positions)),
      all_positions
    )]
    grid_col <- rep("grey70", length(sector_order))
    names(grid_col) <- sector_order
  }

  list(order = sector_order, grid_col = grid_col)
}

#' Named color palette for tRNA structural regions.
#' @noRd
.region_colors <- function() {
  c(
    "acceptor-stem" = "#E41A1C",
    "D-stem" = "#377EB8",
    "D-loop" = "#4DAF4A",
    "AC-stem" = "#984EA3",
    "AC-loop" = "#FF7F00",
    "variable-loop" = "#A65628",
    "T-stem" = "#F781BF",
    "T-loop" = "#999999"
  )
}
