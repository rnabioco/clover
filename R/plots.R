# Plotting functions ----------------------------------------------------------

theme_markdown_axes <- function() {
  rlang::check_installed("ggtext", reason = "to render formatted axis labels.")
  theme(
    axis.title.x = ggtext::element_markdown(),
    axis.title.y = ggtext::element_markdown()
  )
}

# Internal helpers for heatmap -------------------------------------------------

#' Choose text color for contrast against a diverging fill.
#' @noRd
compute_text_color <- function(values, color_limits) {
  threshold <- 0.4 * (color_limits[2] - color_limits[1]) / 2
  ifelse(
    is.na(values) | abs(values) <= threshold,
    "black",
    "white"
  )
}

#' Cluster refs using Ward's D2 on a wide value matrix.
#' @noRd
cluster_refs <- function(data, ref_col, value_col, threshold = NULL) {
  cluster_data <- data

  if (!is.null(threshold)) {
    # Find positions where any row exceeds the threshold
    informative <- cluster_data |>
      dplyr::group_by(sprinzl_label) |>
      dplyr::filter(any(abs(.data[[value_col]]) > threshold)) |>
      dplyr::ungroup()

    if (nrow(informative) > 0) {
      cluster_data <- informative
    }
  }

  wide <- cluster_data |>
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
  rownames(mat)[hc$order]
}

#' Cluster refs within groups, returning ordered refs and group sizes.
#' @return A list with `ref_order` (character) and `group_sizes` (named
#'   integer vector with cumulative counts at each group boundary).
#' @noRd
cluster_refs_by_group <- function(
  data,
  ref_col,
  value_col,
  group_col,
  threshold = NULL
) {
  groups <- unique(data[[group_col]])
  groups <- sort(groups)

  ref_order <- character(0)
  group_sizes <- integer(0)

  for (g in groups) {
    group_data <- data[data[[group_col]] == g, , drop = FALSE]
    refs <- unique(group_data[[ref_col]])
    if (length(refs) > 1) {
      ordered <- cluster_refs(group_data, ref_col, value_col, threshold)
    } else {
      ordered <- refs
    }
    ref_order <- c(ref_order, ordered)
    group_sizes <- c(group_sizes, stats::setNames(length(ref_order), g))
  }

  list(ref_order = ref_order, group_sizes = group_sizes)
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
#' @param label_col Column name (string) with text labels to overlay on
#'   tiles (e.g., nucleotide letters). Default `NULL` (no labels).
#' @param label_min Minimum `abs(value)` to show a label. Default `0.05`.
#' @param label_size Font size for tile labels. Default `2.5`.
#' @param highlight_col Column name (string) of a logical column; `TRUE`
#'   cells get a dot overlay. Default `NULL` (no dots).
#' @param highlight_size Dot size for highlighted cells. Default `0.8`.
#' @param highlight_offset Numeric vector of length 2 giving x/y offsets
#'   from tile center for highlight dots. Default `c(-0.35, 0.35)`.
#' @param cluster_threshold Numeric threshold for noise filtering during
#'   clustering. When non-NULL, only positions where any row has
#'   `abs(value) > cluster_threshold` are used to build the distance
#'   matrix. Falls back to all positions if nothing passes. Default
#'   `NULL`.
#' @param group_col Column name (string) for group-aware clustering. When
#'   provided, rows are clustered within each group and horizontal divider
#'   lines separate groups. Default `NULL`.
#' @param divider_linewidth Line width for group dividers. Default `0.8`.
#' @param fill_name Legend title for the fill scale. Default
#'   `waiver()` (ggplot2 default).
#' @param fill_breaks Numeric vector of legend breaks for the fill
#'   scale. Default `waiver()` (ggplot2 default).
#' @param caption Explanatory text displayed below the plot. Default
#'   `NULL`.
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
  square = TRUE,
  label_col = NULL,
  label_min = 0.05,
  label_size = 2.5,
  highlight_col = NULL,
  highlight_size = 0.4,
  highlight_offset = c(-0.35, 0.35),
  cluster_threshold = NULL,
  group_col = NULL,
  divider_linewidth = 0.8,
  fill_name = waiver(),
  fill_breaks = waiver(),
  caption = NULL
) {
  # --- order x-axis by Sprinzl position ---
  if (!is.factor(data$sprinzl_label)) {
    data$sprinzl_label <- order_sprinzl_positions(data$sprinzl_label)
  }

  # --- cluster rows ---
  refs <- unique(data[[ref_col]])

  if (cluster && length(refs) > 1 && !is.null(group_col)) {
    grouped <- cluster_refs_by_group(
      data,
      ref_col,
      value_col,
      group_col,
      cluster_threshold
    )
    ref_order <- grouped$ref_order
    group_sizes <- grouped$group_sizes
  } else if (cluster && length(refs) > 1) {
    ref_order <- cluster_refs(data, ref_col, value_col, cluster_threshold)
    group_sizes <- NULL
  } else {
    ref_order <- refs
    group_sizes <- NULL
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
      oob = scales::squish,
      name = fill_name,
      breaks = fill_breaks
    ) +
    labs(x = "Sprinzl Position", y = "", caption = caption) +
    cowplot::theme_cowplot() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 5),
      axis.text.y = element_text(size = rel(0.6)),
      legend.position = "bottom",
      legend.justification = c(1, 0),
      legend.direction = "horizontal",
      legend.key.width = grid::unit(1, "cm"),
      legend.key.height = grid::unit(0.3, "cm"),
      legend.text = element_text(size = 7),
      legend.title = element_text(size = 8)
    )

  # --- text labels ---
  if (!is.null(label_col)) {
    plot_data <- dplyr::mutate(
      plot_data,
      .label_display = dplyr::if_else(
        is.na(.data[[value_col]]) | abs(.data[[value_col]]) < label_min,
        NA_character_,
        as.character(.data[[label_col]])
      ),
      .text_color = compute_text_color(.data[[value_col]], color_limits)
    )

    p <- p +
      geom_text(
        data = plot_data,
        aes(label = .data$.label_display, color = .data$.text_color),
        size = label_size,
        fontface = "bold",
        na.rm = TRUE,
        inherit.aes = TRUE,
        show.legend = FALSE
      ) +
      scale_color_identity()
  }

  # --- highlight dots ---
  if (!is.null(highlight_col)) {
    highlight_data <- plot_data[
      !is.na(plot_data[[highlight_col]]) &
        plot_data[[highlight_col]] == TRUE &
        !is.na(plot_data[[value_col]]),
    ]

    if (nrow(highlight_data) > 0) {
      p <- p +
        geom_point(
          data = highlight_data,
          aes(
            x = as.numeric(sprinzl_label) + highlight_offset[1],
            y = as.numeric(.data[[ref_col]]) + highlight_offset[2]
          ),
          color = "black",
          size = highlight_size,
          inherit.aes = FALSE
        )
    }
  }

  # --- group dividers ---
  if (!is.null(group_sizes)) {
    n_refs <- length(ref_order)
    for (i in seq_along(group_sizes)[-length(group_sizes)]) {
      # y-axis is reversed, so divider is at n_refs - boundary + 0.5
      boundary <- group_sizes[i]
      p <- p +
        geom_hline(
          yintercept = n_refs - boundary + 0.5,
          color = "black",
          linewidth = divider_linewidth
        )
    }
  }

  # --- caption theme ---
  if (!is.null(caption)) {
    p <- p +
      theme(
        plot.caption = element_text(
          hjust = 0,
          size = 10,
          margin = margin(t = 8)
        ),
        plot.caption.position = "plot"
      )
  }

  if (square) {
    p <- p + coord_fixed(ratio = 1)
  }

  p
}

#' Plot a volcano plot of differential expression results.
#'
#' Creates a volcano plot from the tibble returned by
#' [tidy_deseq_results()]. Significant points are labeled with
#' [ggrepel::geom_text_repel()].
#'
#' @param data A tibble from [tidy_deseq_results()] with at least
#'   `log2FoldChange`, `padj`, and `significant` columns.
#' @param lab_col Column name (string) used for point labels. Default
#'   `"ref"`.
#' @param padj_cutoff Numeric; draws a dashed horizontal line at
#'   `-log10(padj_cutoff)`. Default `0.05`.
#' @param max_overlaps Maximum number of overlapping labels passed to
#'   [ggrepel::geom_text_repel()]. Default `20`.
#' @param point_size Numeric size for [ggplot2::geom_point()]. Default
#'   `1.5`.
#' @param label_size Numeric size for [ggrepel::geom_text_repel()].
#'   Default `3`.
#' @param sig_color Color for significant points. Default `"#D55E00"`.
#' @param nonsig_color Color for non-significant points. Default
#'   `"grey60"`.
#'
#' @return A ggplot object.
#'
#' @export
#'
#' @examples
#' res <- tibble::tibble(
#'   ref = paste0("tRNA-", 1:10),
#'   log2FoldChange = rnorm(10),
#'   pvalue = c(rep(0.001, 3), rep(0.5, 7)),
#'   padj = c(rep(0.01, 3), rep(0.8, 7)),
#'   significant = c(rep(TRUE, 3), rep(FALSE, 7))
#' )
#' plot_volcano(res)
plot_volcano <- function(
  data,
  lab_col = "ref",
  padj_cutoff = 0.05,
  max_overlaps = 20,
  point_size = 1.5,
  label_size = 3,
  sig_color = "#D55E00",
  nonsig_color = "grey60"
) {
  rlang::check_installed("ggrepel", reason = "to label significant points.")

  p <- ggplot(data, aes(x = log2FoldChange, y = -log10(padj))) +
    geom_point(
      aes(color = significant),
      size = point_size,
      alpha = 0.7
    ) +
    ggrepel::geom_text_repel(
      data = function(x) dplyr::filter(x, significant),
      aes(label = .data[[lab_col]]),
      size = label_size,
      max.overlaps = max_overlaps
    ) +
    scale_color_manual(
      values = stats::setNames(
        c(nonsig_color, sig_color),
        c(FALSE, TRUE)
      )
    ) +
    geom_vline(xintercept = 0, color = "grey80") +
    geom_hline(
      yintercept = -log10(padj_cutoff),
      linetype = "dashed",
      color = "grey40"
    ) +
    labs(
      x = "log<sub>2</sub> Fold Change",
      y = "\u2212log<sub>10</sub>(p-value)"
    ) +
    cowplot::theme_cowplot() +
    theme(legend.position = "none") +
    theme_markdown_axes()

  p
}

#' Plot abundance changes versus charging ratio changes.
#'
#' Creates a scatter plot comparing tRNA abundance changes (from DESeq2)
#' with charging ratio changes on a single plot. Significant points are
#' colored by quadrant and labeled with
#' [ggrepel::geom_text_repel()].
#'
#' @param deseq_res A tibble from [tidy_deseq_results()] with at least
#'   `ref`, `log2FoldChange`, and `padj` columns.
#' @param charging_diffs A tibble from [compute_charging_diffs()] with
#'   at least `ref` and `diff` columns.
#' @param lab_col Column name (string) used for point labels. Default
#'   `"ref"`.
#' @param padj_cutoff Numeric; significance threshold for `padj`.
#'   Default `0.05`.
#' @param max_overlaps Maximum number of overlapping labels passed to
#'   [ggrepel::geom_text_repel()]. Default `20`.
#' @param point_size Numeric size for [ggplot2::geom_point()]. Default
#'   `2`.
#' @param label_size Numeric size for [ggrepel::geom_text_repel()].
#'   Default `3`.
#'
#' @return A ggplot object.
#'
#' @export
#'
#' @examples
#' deseq_res <- tibble::tibble(
#'   ref = paste0("tRNA-", 1:6),
#'   log2FoldChange = c(1, -1, 0.5, -0.5, 2, -2),
#'   padj = c(0.01, 0.02, 0.5, 0.6, 0.001, 0.003)
#' )
#' charging_diffs <- tibble::tibble(
#'   ref = paste0("tRNA-", 1:6),
#'   diff = c(0.1, -0.1, 0.05, -0.05, -0.2, 0.15),
#'   se_diff = rep(0.03, 6)
#' )
#' plot_abundance_charging(deseq_res, charging_diffs)
plot_abundance_charging <- function(
  deseq_res,
  charging_diffs,
  lab_col = "ref",
  padj_cutoff = 0.05,
  max_overlaps = 20,
  point_size = 2,
  label_size = 3
) {
  rlang::check_installed("ggrepel", reason = "to label significant points.")

  data <- dplyr::inner_join(deseq_res, charging_diffs, by = "ref")

  data <- dplyr::mutate(
    data,
    significant = !is.na(.data$padj) & .data$padj < padj_cutoff,
    quadrant = dplyr::case_when(
      !significant ~ "ns",
      log2FoldChange >= 0 & diff >= 0 ~ "up_up",
      log2FoldChange < 0 & diff < 0 ~ "down_down",
      log2FoldChange >= 0 & diff < 0 ~ "up_down",
      log2FoldChange < 0 & diff >= 0 ~ "down_up"
    )
  )

  quad_colors <- c(
    up_up = "#D55E00",
    down_down = "#0072B2",
    up_down = "#CC79A7",
    down_up = "#009E73",
    ns = "grey60"
  )

  quad_labels <- c(
    up_up = "Up / Up",
    down_down = "Down / Down",
    up_down = "Up / Down",
    down_up = "Down / Up",
    ns = "NS"
  )

  # Only include quadrants present in the data
  present <- intersect(names(quad_colors), unique(data$quadrant))

  p <- ggplot(data, aes(x = log2FoldChange, y = diff)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey40") +
    geom_point(
      aes(color = quadrant),
      size = point_size,
      alpha = 0.7
    ) +
    ggrepel::geom_text_repel(
      data = function(x) dplyr::filter(x, .data$significant),
      aes(label = .data[[lab_col]]),
      size = label_size,
      max.overlaps = max_overlaps
    ) +
    scale_color_manual(
      values = quad_colors[present],
      labels = quad_labels[present]
    ) +
    labs(
      x = "log<sub>2</sub> Fold Change (abundance)",
      y = "Charging ratio difference",
      color = "Abundance / Charging"
    ) +
    cowplot::theme_cowplot() +
    theme(legend.position = "bottom") +
    theme_markdown_axes()

  p
}

#' Plot per-tRNA charging ratio differences.
#'
#' Create a dot plot with error bars showing the difference in charging
#' ratio between two conditions for each tRNA. Consumes the tibble
#' returned by [compute_charging_diffs()].
#'
#' @param data A tibble from [compute_charging_diffs()] with at least
#'   `ref` (factor), `diff`, and `se_diff` columns.
#' @param point_size Numeric size for [ggplot2::geom_point()]. Default
#'   `2.5`.
#'
#' @return A ggplot object.
#'
#' @export
#'
#' @examples
#' df <- tibble::tibble(
#'   ref = forcats::fct_inorder(paste0("tRNA-", 1:5)),
#'   diff = c(-0.1, -0.05, 0.02, 0.08, 0.15),
#'   se_diff = rep(0.03, 5)
#' )
#' plot_charging_diffs(df)
plot_charging_diffs <- function(data, point_size = 2.5) {
  ggplot(data, aes(x = diff, y = ref)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "gray50") +
    geom_point(size = point_size) +
    geom_linerange(aes(xmin = diff - se_diff, xmax = diff + se_diff)) +
    labs(
      x = "Difference in charging ratio",
      y = ""
    ) +
    cowplot::theme_minimal_vgrid()
}

#' Plot per-position base-calling error profiles.
#'
#' Create a line plot of per-position base-calling error rates, faceted
#' by tRNA and colored by condition. Optionally overlay vertical dashed
#' lines at known modification positions.
#'
#' @param data A summarized bcerror tibble with columns `ref`, `pos`,
#'   `condition`, and `mean_error`.
#' @param refs Character vector of tRNA names to plot (filters
#'   `data$ref`). If `NULL` (default), all tRNAs are plotted.
#' @param mods Optional tibble with `ref` and `pos` columns marking
#'   modification positions (e.g., from [fetch_modomics_mods()]).
#' @param colors Named character vector of colors for conditions.
#'   Default `c(ctl = "#0072B2", inf = "#D55E00")`.
#' @param ncol Number of columns for [ggplot2::facet_wrap()]. Default
#'   `1`.
#'
#' @return A ggplot object.
#'
#' @export
#'
#' @examples
#' df <- tidyr::expand_grid(
#'   ref = c("tRNA-Ala", "tRNA-Gly"),
#'   pos = 1:20,
#'   condition = c("ctl", "inf")
#' )
#' df$mean_error <- runif(nrow(df), 0, 0.3)
#' plot_bcerror_profile(df)
plot_bcerror_profile <- function(
  data,
  refs = NULL,
  mods = NULL,
  colors = c(ctl = "#0072B2", inf = "#D55E00"),
  ncol = 1
) {
  if (!is.null(refs)) {
    data <- dplyr::filter(data, ref %in% refs)
  }

  p <- ggplot(data, aes(x = pos, y = mean_error, color = condition)) +
    geom_line(linewidth = 0.5) +
    geom_point(size = 0.8)

  if (!is.null(mods)) {
    if (!is.null(refs)) {
      mods <- dplyr::filter(mods, ref %in% refs)
    }
    p <- p +
      geom_vline(
        data = mods,
        aes(xintercept = pos),
        linetype = "dashed",
        color = "grey40",
        alpha = 0.5,
        inherit.aes = FALSE
      )
  }

  p +
    facet_wrap(~ref, ncol = ncol, scales = "free_y") +
    scale_color_manual(values = colors) +
    labs(
      x = "Position",
      y = "Mean base-calling error rate"
    ) +
    cowplot::theme_minimal_hgrid() +
    theme(legend.position = "top")
}

# Modification landscape helpers -----------------------------------------------

#' Create region shading layers from position-region data.
#' @return A list of `geom_rect()` layers.
#' @noRd
add_region_shading <- function(data, pos_col, region_col) {
  region_data <- data |>
    dplyr::filter(!is.na(.data[[region_col]])) |>
    dplyr::distinct(.data[[pos_col]], .data[[region_col]]) |>
    dplyr::arrange(.data[[pos_col]])

  if (nrow(region_data) == 0) {
    return(list())
  }

  region_data <- dplyr::mutate(
    region_data,
    .region_change = .data[[region_col]] !=
      dplyr::lag(
        .data[[region_col]],
        default = ""
      ),
    .seg_id = cumsum(.data$.region_change)
  )

  region_segs <- region_data |>
    dplyr::group_by(.data$.seg_id, .data[[region_col]]) |>
    dplyr::summarise(
      .xmin = min(.data[[pos_col]]) - 0.5,
      .xmax = max(.data[[pos_col]]) + 0.5,
      .groups = "drop"
    )

  palette <- region_colors()

  purrr::map(seq_len(nrow(region_segs)), function(i) {
    seg <- region_segs[i, ]
    fill <- unname(
      palette[match(seg[[region_col]], names(palette))]
    )
    if (is.na(fill)) {
      fill <- "grey70"
    }
    geom_rect(
      data = seg,
      aes(xmin = .data$.xmin, xmax = .data$.xmax, ymin = -Inf, ymax = Inf),
      fill = fill,
      alpha = 0.15,
      inherit.aes = FALSE
    )
  })
}

#' Build a secondary x-axis with Sprinzl labels.
#' @return A list with `$scale` and `$theme` elements to add to a ggplot.
#' @noRd
create_sprinzl_axis <- function(data, pos_col, sprinzl_col, mod_col = NULL) {
  pos_map <- data |>
    dplyr::distinct(.data[[pos_col]], .data[[sprinzl_col]]) |>
    dplyr::filter(!is.na(.data[[sprinzl_col]])) |>
    dplyr::arrange(.data[[pos_col]])

  sec_breaks <- pos_map[[pos_col]]
  sec_labels <- as.character(pos_map[[sprinzl_col]])

  have_ggtext <- rlang::is_installed("ggtext")

  if (!is.null(mod_col) && have_ggtext) {
    mod_map <- data |>
      dplyr::distinct(.data[[pos_col]], .data[[mod_col]]) |>
      dplyr::filter(!is.na(.data[[mod_col]]))

    mod_positions <- mod_map[[pos_col]][
      !is.na(mod_map[[mod_col]]) & mod_map[[mod_col]] == TRUE
    ]

    sec_labels <- ifelse(
      sec_breaks %in% mod_positions,
      sprintf("<span style='color:#D55E00'>**%s**</span>", sec_labels),
      sec_labels
    )

    theme_el <- ggtext::element_markdown(
      angle = 45,
      size = 7,
      hjust = 0
    )
  } else {
    theme_el <- element_text(angle = 45, size = 7, hjust = 0)
  }

  list(
    scale = scale_x_continuous(
      sec.axis = ggplot2::sec_axis(
        ~.,
        breaks = sec_breaks,
        labels = sec_labels,
        name = "Sprinzl Position"
      )
    ),
    theme = theme(axis.text.x.top = theme_el)
  )
}

#' Plot per-tRNA modification landscape profiles.
#'
#' Create a stacked panel plot showing multiple metrics along the tRNA
#' sequence. Each metric gets its own panel sharing a common x-axis.
#' Optionally adds structural region background shading and a secondary
#' x-axis with Sprinzl position labels.
#'
#' @param data A data frame with a position column and one or more
#'   metric columns to plot.
#' @param metrics Character vector of column names to plot as stacked
#'   panels (one panel per metric).
#' @param pos_col Column name (string) for x-axis positions. Default
#'   `"pos"`.
#' @param region_col Optional column name (string) for structural
#'   region labels, used for background shading. Default `NULL`.
#' @param sprinzl_col Optional column name (string) for Sprinzl
#'   position labels shown on a secondary x-axis. Default `NULL`.
#' @param mod_col Optional column name (string) of a logical column;
#'   `TRUE` positions are highlighted on the Sprinzl axis (requires
#'   ggtext). Default `NULL`.
#' @param title Plot title. Default `NULL`.
#' @param heights Numeric vector of relative panel heights. Default
#'   `NULL` (equal heights).
#'
#' @return A patchwork object combining stacked ggplot panels.
#'
#' @export
#'
#' @examples
#' df <- data.frame(
#'   pos = rep(1:20, 2),
#'   condition = rep(c("ctl", "mut"), each = 20),
#'   error_rate = runif(40, 0, 0.3),
#'   signal = rnorm(40, sd = 0.1)
#' )
#' plot_mod_landscape(df, metrics = c("error_rate", "signal"))
plot_mod_landscape <- function(
  data,
  metrics,
  pos_col = "pos",
  region_col = NULL,
  sprinzl_col = NULL,
  mod_col = NULL,
  title = NULL,
  heights = NULL
) {
  rlang::check_installed(
    "patchwork",
    reason = "to stack landscape panels."
  )

  region_layers <- if (!is.null(region_col)) {
    add_region_shading(data, pos_col, region_col)
  } else {
    list()
  }

  n_metrics <- length(metrics)
  panels <- vector("list", n_metrics)

  for (i in seq_along(metrics)) {
    metric <- metrics[i]
    is_first <- i == 1
    is_last <- i == n_metrics

    p <- ggplot(data, aes(x = .data[[pos_col]], y = .data[[metric]])) +
      region_layers +
      geom_line(linewidth = 0.5) +
      geom_point(size = 0.8) +
      labs(y = metric) +
      cowplot::theme_minimal_hgrid()

    if (!is_last) {
      p <- p +
        labs(x = NULL) +
        theme(axis.text.x = element_blank())
    } else {
      p <- p + labs(x = "Position")
    }

    # Add Sprinzl secondary axis to top panel
    if (is_first && !is.null(sprinzl_col)) {
      sprinzl_axis <- create_sprinzl_axis(
        data,
        pos_col,
        sprinzl_col,
        mod_col
      )
      p <- p + sprinzl_axis$scale + sprinzl_axis$theme
    }

    panels[[i]] <- p
  }

  combined <- patchwork::wrap_plots(panels, ncol = 1, heights = heights)

  if (!is.null(title)) {
    combined <- combined +
      patchwork::plot_annotation(title = title)
  }

  combined
}

#' Plot PCoA of tRNA rewiring scores.
#'
#' Create a scatter plot of PCoA coordinates from
#' [perform_pcoa()], with points sized by the number of non-zero
#' comparisons and colored by Euclidean rewiring magnitude. The top
#' `n_label` isodecoders are labeled with
#' [ggrepel::geom_text_repel()].
#'
#' @param pcoa_result A list from [perform_pcoa()] with elements
#'   `coordinates` and `variance_explained`.
#' @param rewiring_scores A tibble from [calculate_rewiring_scores()]
#'   with at least `isodecoder`, `euclidean_magnitude`, and
#'   `n_nonzero` columns.
#' @param title Plot title. Default `"PCoA of tRNA rewiring"`.
#' @param n_label Number of top isodecoders to label. Default `10`.
#'
#' @return A ggplot object.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' mat <- prepare_rewiring_matrix(ror_data)
#' scores <- calculate_rewiring_scores(mat)
#' pcoa <- perform_pcoa(mat)
#' plot_pcoa_rewiring(pcoa, scores)
#' }
plot_pcoa_rewiring <- function(
  pcoa_result,
  rewiring_scores,
  title = "PCoA of tRNA rewiring",
  n_label = 10
) {
  rlang::check_installed("ggrepel", reason = "to label top isodecoders.")

  plot_data <- dplyr::left_join(
    pcoa_result$coordinates,
    rewiring_scores,
    by = "isodecoder"
  )

  top_trnas <- plot_data |>
    dplyr::arrange(dplyr::desc(euclidean_magnitude)) |>
    utils::head(n_label) |>
    dplyr::pull(isodecoder)

  plot_data <- dplyr::mutate(
    plot_data,
    label = ifelse(isodecoder %in% top_trnas, isodecoder, NA_character_)
  )

  ggplot(plot_data, aes(x = PC1, y = PC2)) +
    geom_point(
      aes(size = n_nonzero, color = euclidean_magnitude),
      alpha = 0.7
    ) +
    ggrepel::geom_text_repel(
      aes(label = label),
      size = 3,
      max.overlaps = 20,
      box.padding = 0.5,
      segment.color = "grey50"
    ) +
    scale_color_viridis_c(option = "plasma", name = "Rewiring\nmagnitude") +
    scale_size_continuous(name = "Non-zero\nchanges", range = c(2, 8)) +
    labs(
      title = title,
      x = paste0(
        "PC1 (",
        round(pcoa_result$variance_explained[1], 1),
        "%)"
      ),
      y = paste0(
        "PC2 (",
        round(pcoa_result$variance_explained[2], 1),
        "%)"
      )
    ) +
    cowplot::theme_cowplot() +
    theme(legend.position = "right")
}
