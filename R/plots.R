# Plotting functions ----------------------------------------------------------

#' Plot base-calling error in a heatmap.
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

#' Plot a delta-signal modification heatmap.
#'
#' Create a diverging heatmap of modification signal changes (e.g., mutant
#' minus wild-type) across tRNA families and Sprinzl positions. Rows can be
#' optionally clustered using Ward's D2 hierarchical clustering.
#'
#' @param data A data frame with at least three columns: one for tRNA
#'   family/reference (y-axis), one for Sprinzl position labels (x-axis),
#'   and one for the fill value.
#' @param value_col Column name (string) for fill values. Default `"value"`.
#' @param ref_col Column name (string) for tRNA families (y-axis). Default
#'   `"ref"`.
#' @param cluster Logical; cluster rows with Ward's D2? Default `TRUE`.
#' @param color_limits Numeric vector of length 2 giving symmetric limits
#'   for the color scale. Default `c(-0.25, 0.25)`.
#' @param color_low Color for negative values. Default `"#0072B2"` (blue).
#' @param color_high Color for positive values. Default `"#D55E00"` (red).
#' @param na_value Color for missing positions. Default `"gray80"`.
#' @param square Logical; use `coord_fixed(ratio = 1)`? Default `TRUE`.
#'
#' @return A ggplot object.
#'
#' @export
#'
#' @examples
#' df <- tidyr::expand_grid(
#'   ref = paste0("tRNA-", c("Ala", "Gly", "Ser")),
#'   sprinzl_label = as.character(1:10)
#' )
#' df$value <- rnorm(nrow(df), sd = 0.1)
#' plot_mod_heatmap(df)
plot_mod_heatmap <- function(
  data,
  value_col = "value",
  ref_col = "ref",
  cluster = TRUE,
  color_limits = c(-0.25, 0.25),
  color_low = "#0072B2",
  color_high = "#D55E00",
  na_value = "gray80",
  square = TRUE
) {
  # --- order x-axis by Sprinzl position ---
  data$sprinzl_label <- order_sprinzl_positions(data$sprinzl_label)

  # --- cluster rows ---
  refs <- unique(data[[ref_col]])

  if (cluster && length(refs) > 1) {
    wide <- data |>
      dplyr::select(
        dplyr::all_of(c(ref_col, "sprinzl_label", value_col))
      ) |>
      tidyr::pivot_wider(
        names_from = sprinzl_label,
        values_from = dplyr::all_of(value_col),
        values_fill = 0
      ) |>
      as.data.frame()

    rownames(wide) <- wide[[ref_col]]
    mat <- as.matrix(wide[, -1, drop = FALSE])
    mat[is.na(mat)] <- 0

    hc <- stats::hclust(stats::dist(mat), method = "ward.D2")
    ref_order <- rownames(mat)[hc$order]
  } else {
    ref_order <- refs
  }

  # --- complete grid so missing cells show as gray ---
  all_positions <- levels(data$sprinzl_label)
  complete_grid <- tidyr::expand_grid(
    !!ref_col := ref_order,
    sprinzl_label = all_positions
  )

  plot_data <- dplyr::left_join(
    complete_grid,
    data,
    by = c(ref_col, "sprinzl_label")
  ) |>
    dplyr::mutate(
      sprinzl_label = factor(sprinzl_label, levels = all_positions),
      !!ref_col := factor(.data[[ref_col]], levels = rev(ref_order))
    )

  # --- build plot ---
  p <- ggplot(
    plot_data,
    aes(
      x = sprinzl_label,
      y = .data[[ref_col]],
      fill = .data[[value_col]]
    )
  ) +
    geom_tile(color = "white", linewidth = 0.3) +
    scale_fill_gradient2(
      low = color_low,
      mid = "white",
      high = color_high,
      midpoint = 0,
      na.value = na_value,
      limits = color_limits,
      oob = scales::squish
    ) +
    labs(x = "Sprinzl Position", y = "") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
      panel.grid = element_blank(),
      legend.position = "bottom",
      legend.key.width = grid::unit(1.5, "cm")
    )

  if (square) {
    p <- p + coord_fixed(ratio = 1)
  }

  p
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
