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
