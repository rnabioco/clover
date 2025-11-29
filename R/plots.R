# Plotting functions ----------------------------------------------------------

#' Internal helper to create heatmap plot
#'
#' @param plot_data Filtered bcerror data
#' @param value Column name for fill
#' @param show_regions Logical for region annotations
#' @param title_suffix Optional title suffix (e.g., "Type I", "Type II")
#'
#' @noRd
.plot_heatmap_internal <- function(
    plot_data,
    value,
    show_regions,
    title_suffix = NULL
) {
  # Filter to positions where at least one tRNA has non-NA data
  positions_with_data <- plot_data |>
    dplyr::filter(!is.na(.data[[value]])) |>
    dplyr::pull(global_index) |>
    unique()

  plot_data <- plot_data |>
    dplyr::filter(global_index %in% positions_with_data)

  # Get axis labels for positions with data
  all_indices <- sort(unique(plot_data$global_index))

  labels_lookup <- plot_data |>
    dplyr::distinct(global_index, sprinzl_label) |>
    dplyr::mutate(
      sprinzl_label = dplyr::if_else(
        is.na(sprinzl_label) | sprinzl_label == "-1",
        "",
        sprinzl_label
      )
    )

  # Create label vector for positions with data
  axis_labels <- sapply(all_indices, function(idx) {
    label <- labels_lookup$sprinzl_label[labels_lookup$global_index == idx]
    if (length(label) == 0) "" else label[1]
  })

  # Compute region bounds BEFORE converting global_index to factor
  if (show_regions) {
    regions <- get_region_bounds(plot_data)

    # Convert numeric positions to factor-based positions for discrete scale
    # Factor levels are at positions 1, 2, 3, ... so we map global_index to those
    index_to_position <- setNames(seq_along(all_indices), all_indices)

    regions <- regions |>
      dplyr::filter(start %in% all_indices, end %in% all_indices) |>
      dplyr::mutate(
        xmin = index_to_position[as.character(start)] - 0.5,
        xmax = index_to_position[as.character(end)] + 0.5
      )
  }

  # Convert global_index to factor to use discrete scale (omits empty positions)
  plot_data <- plot_data |>
    dplyr::mutate(global_index = factor(global_index, levels = all_indices))

  # Simplify tRNA labels: strip "nuc-tRNA-" prefix if no mito tRNAs present
  has_mito <- any(grepl("^mito-", plot_data$ref))
  if (!has_mito) {
    plot_data <- plot_data |>
      dplyr::mutate(ref = sub("^nuc-tRNA-", "", ref))
  }

  # Build plot
  p <- ggplot(
    plot_data,
    aes(
      x = global_index,
      y = forcats::fct_rev(ref),
      fill = .data[[value]]
    )
  ) +
    geom_tile() +
    coord_fixed(ratio = 1) +
    scale_fill_viridis_c(
      option = "magma",
      na.value = "grey90",
      name = value
    ) +
    scale_x_discrete(
      labels = axis_labels,
      expand = c(0, 0)
    ) +
    labs(
      x = "Position (Sprinzl)",
      y = NULL
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 7),
      axis.text.y = element_text(size = 6),
      panel.grid = element_blank(),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA)
    )

  # Add title if suffix provided
  if (!is.null(title_suffix)) {
    p <- p + labs(title = paste("tRNAs -", title_suffix))
  }

  # Add region annotations if requested
  if (show_regions && nrow(regions) > 0) {
    p <- p +
      geom_rect(
        data = regions,
        aes(
          xmin = xmin,
          xmax = xmax,
          ymin = -Inf,
          ymax = Inf,
          fill = NULL
        ),
        alpha = 0,
        color = "grey50",
        linewidth = 0.25,
        inherit.aes = FALSE
      )
  }

  p
}


#' Plot base-calling error heatmap with global coordinates
#'
#' Creates a heatmap visualization of base-calling error rates across tRNA
#' positions. When data contains offset and type columns (from
#' [add_global_coords()]), splits into separate panels for each coordinate
#' group to ensure proper Sprinzl position alignment.
#'
#' @param bcerror Tibble with base-calling error data. Must have global
#'   coordinates added via [add_global_coords()].
#' @param value Column name to map to heatmap color. Default is "error_rate".
#'   Can use "mis", "ins", "del", or any numeric column.
#' @param show_regions Logical. If TRUE, adds grey boxes around structural
#'   regions (acceptor-stem, D-loop, etc.).
#' @param split_by_group Logical. If TRUE (default), returns stacked plots for
#'   each offset×type group. If FALSE, returns single plot with all tRNAs.
#'
#' @return A patchwork object with separate heatmaps for each coordinate group
#'   stacked vertically (if split_by_group=TRUE), or a single ggplot object
#'   (if FALSE).
#'
#'   All Sprinzl position labels are shown on the x-axis. Positions without
#'   Sprinzl labels display "NA".
#'
#' @export
#'
#' @examples
#' bcerr <- read_bcerror(clover_example("yeast/grande.bcerr.tsv.gz")) |>
#'   add_global_coords("sacCer")
#'
#' # Get stacked plots for each offset×type group (default)
#' plot_bcerror_heatmap(bcerr)
#'
#' # Single plot with all tRNAs (may have alignment issues across groups)
#' plot_bcerror_heatmap(bcerr, split_by_group = FALSE)
plot_bcerror_heatmap <- function(
    bcerror,
    value = "error_rate",
    show_regions = TRUE,
    split_by_group = TRUE
) {
  if (!"global_index" %in% names(bcerror)) {
    stop("bcerror must have global coordinates. Use add_global_coords() first.")
  }

  # Filter to rows with valid global coordinates
  plot_data <- bcerror |>
    dplyr::filter(!is.na(global_index))

  if (split_by_group && "offset" %in% names(plot_data) && "type" %in% names(plot_data)) {
    # Split by offset×type groups for proper alignment
    groups <- plot_data |>
      dplyr::distinct(offset, type) |>
      dplyr::arrange(type, offset)

    # Create a plot for each group
    plots <- lapply(seq_len(nrow(groups)), function(i) {
      grp_offset <- groups$offset[i]
      grp_type <- groups$type[i]

      grp_data <- plot_data |>
        dplyr::filter(offset == grp_offset, type == grp_type)

      n_trnas <- length(unique(grp_data$ref))

      # Create descriptive title
      type_label <- ifelse(grp_type == "type1", "Type I", "Type II")
      title <- sprintf("%s (offset %s) - %d tRNAs",
                       type_label,
                       ifelse(grp_offset >= 0, paste0("+", grp_offset), grp_offset),
                       n_trnas)

      .plot_heatmap_internal(grp_data, value, show_regions, title)
    })

    # Stack vertically using patchwork
    patchwork::wrap_plots(plots, ncol = 1)
  } else if (split_by_group) {
    # Fallback to Type I/II split if offset/type columns missing
    plot_data <- plot_data |>
      dplyr::mutate(trna_type = classify_trna_type(ref))

    type1_data <- plot_data |> dplyr::filter(trna_type == "Type I")
    type2_data <- plot_data |> dplyr::filter(trna_type == "Type II")

    p_type1 <- .plot_heatmap_internal(type1_data, value, show_regions, "Type I")
    p_type2 <- .plot_heatmap_internal(type2_data, value, show_regions, "Type II")

    patchwork::wrap_plots(p_type1, p_type2, ncol = 1)
  } else {
    # Single plot (old behavior)
    .plot_heatmap_internal(plot_data, value, show_regions, NULL)
  }
}

#' Plot base-calling error in a heatmap (legacy)
#'
#' @param tbl tibble of bc-delta values
#' @param data sequence to structure file
#' @param title_suffix string
#' @param include_legend Logical, include legend in the plot.
#'
#' @export
plot_bcerror <- function(tbl, data, title_suffix, include_legend = TRUE) {
  data <-
    dplyr::mutate(
      data,
      pos = factor(as.numeric(data$pos), levels = 1:96),
      ref = forcats::fct_rev(data$ref)
    )

  complete_grid <- tidyr::expand_grid(
    pos = levels(data$pos),
    ref = levels(data$ref)
  )

  full_data <- dplyr::left_join(
    complete_grid,
    data,
    by = c("pos", "ref")
  )

  ggplot(
    full_data,
    aes(
      x = pos,
      y = ref,
      fill = BC_delta
    )
  ) +
    geom_tile(
      color = "white",
      size = 0.1
    ) +
    scale_fill_gradient2(
      low = "#0072B2",
      high = "#D55E00",
      mid = "white",
      midpoint = 0,
      limits = c(-0.5, 0.5),
      na.value = "#D3D3D3",
      breaks = seq(-0.5, 0.5, by = 0.25)
    ) +
    scale_x_discrete(
      labels = trna_consensus_labels
    ) +
    labs(
      title = paste("Difference in basecalling error", title_suffix),
      x = "Position",
      y = ""
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, size = 8),
      legend.pos = ifelse(include_legend, "bottom", "none")
    )
}

# Constants ---------------------------------------------------------------

#' Labels for consensus tRNA secondary structure
#' @export
trna_consensus_labels <- c(
  0:17,
  "17a",
  18:20,
  "20a",
  "20b",
  21:48,
  # variable loop
  paste0("e", c(11:17, 1:5, 27:21)),
  49:73
)

#' Amino acid properties
#' @export
# fmt: skip
aa_props <- tibble::tribble(
  ~aa_full,        ~aa3,  ~aa1,  ~color,    ~mw,
  "Alanine",       "Ala", "A",   "#8CFF8C", 89.09,
  "Arginine",      "Arg", "R",   "#00007C", 174.20,
  "Asparagine",    "Asn", "N",   "#FF7C70", 132.12,
  "Aspartate",     "Asp", "D",   "#A00042", 133.10,
  "Cysteine",      "Cys", "C",   "#FFFF70", 121.15,
  "Glutamate",     "Glu", "E",   "#A00042", 147.13,
  "Glutamine",     "Gln", "Q",   "#FF4C4C", 146.15,
  "Glycine",       "Gly", "G",   "#FFFFFF", 75.07,
  "Histidine",     "His", "H",   "#7070FF", 155.16,
  "Isoleucine",    "Ile", "I",   "#004C00", 131.17,
  "Leucine",       "Leu", "L",   "#455E45", 131.17,
  "Lysine",        "Lys", "K",   "#4747B8", 146.19,
  "Methionine",    "Met", "M",   "#B8A042", 149.21,
  "Phenylalanine", "Phe", "F",   "#534C52", 165.19,
  "Proline",       "Pro", "P",   "#525252", 115.13,
  "Serine",        "Ser", "S",   "#FF7042", 105.09,
  "Threonine",     "Thr", "T",   "#B84C00", 119.12,
  "Tryptophan",    "Trp", "W",   "#4F4600", 204.23,
  "Tyrosine",      "Tyr", "Y",   "#8C704C", 181.19,
  "Valine",        "Val", "V",   "#FF8CFF", 117.15
)
