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
  # Get ALL axis labels for all global indices
  all_indices <- sort(unique(plot_data$global_index))

  labels_lookup <- plot_data |>
    dplyr::distinct(global_index, sprinzl_label) |>
    dplyr::mutate(
      sprinzl_label = ifelse(sprinzl_label == "-1", "NA", sprinzl_label)
    )

  # Create full label vector with NA for missing positions
  axis_labels <- sapply(all_indices, function(idx) {
    label <- labels_lookup$sprinzl_label[labels_lookup$global_index == idx]
    if (length(label) == 0) "NA" else label[1]
  })
  names(axis_labels) <- all_indices

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
    scale_fill_viridis_c(
      option = "magma",
      na.value = "grey90",
      name = value
    ) +
    scale_x_continuous(
      breaks = as.numeric(names(axis_labels)),
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
  if (show_regions) {
    regions <- get_region_bounds(plot_data)

    p <- p +
      geom_rect(
        data = regions,
        aes(
          xmin = start - 0.5,
          xmax = end + 0.5,
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
#' positions. Splits into separate plots for Type I (standard) and Type II
#' (extended variable loop) tRNAs, stacked vertically.
#'
#' Type II tRNAs (Leu, Ser, Tyr, SeC) have extended variable loops with 9-24
#' extra nucleotides and are displayed in a separate heatmap.
#'
#' @param bcerror Tibble with base-calling error data. Must have global
#'   coordinates added via [add_global_coords()].
#' @param value Column name to map to heatmap color. Default is "error_rate".
#'   Can use "mis", "ins", "del", or any numeric column.
#' @param show_regions Logical. If TRUE, adds grey boxes around structural
#'   regions (acceptor-stem, D-loop, etc.).
#' @param split_by_type Logical. If TRUE (default), returns stacked plots for
#'   Type I and Type II tRNAs. If FALSE, returns single plot with all tRNAs.
#'
#' @return A patchwork object with Type I and Type II heatmaps stacked
#'   vertically (if split_by_type=TRUE), or a single ggplot object (if FALSE).
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
#' # Get stacked plots for Type I and Type II (default)
#' plot_bcerror_heatmap(bcerr)
#'
#' # Single plot with all tRNAs
#' plot_bcerror_heatmap(bcerr, split_by_type = FALSE)
plot_bcerror_heatmap <- function(
    bcerror,
    value = "error_rate",
    show_regions = TRUE,
    split_by_type = TRUE
) {
  if (!"global_index" %in% names(bcerror)) {
    stop("bcerror must have global coordinates. Use add_global_coords() first.")
  }

  # Filter to rows with valid global coordinates
  plot_data <- bcerror |>
    dplyr::filter(!is.na(global_index))

  # Add tRNA type classification
  plot_data <- plot_data |>
    dplyr::mutate(
      trna_type = classify_trna_type(ref)
    )

  if (split_by_type) {
    # Split into Type I and Type II
    type1_data <- plot_data |> dplyr::filter(trna_type == "Type I")
    type2_data <- plot_data |> dplyr::filter(trna_type == "Type II")

    # Create two separate plots
    p_type1 <- .plot_heatmap_internal(type1_data, value, show_regions, "Type I")
    p_type2 <- .plot_heatmap_internal(type2_data, value, show_regions, "Type II")

    # Stack vertically using patchwork
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
