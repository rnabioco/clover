# Plotting functions ----------------------------------------------------------

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
    cowplot::theme_cowplot() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
      legend.position = "bottom",
      legend.key.width = grid::unit(1.5, "cm")
    )

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
#'   `log2FoldChange`, `pvalue`, and `significant` columns.
#' @param lab_col Column name (string) used for point labels. Default
#'   `"tRNA"`.
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
#'   tRNA = paste0("tRNA-", 1:10),
#'   log2FoldChange = rnorm(10),
#'   pvalue = c(rep(0.001, 3), rep(0.5, 7)),
#'   padj = c(rep(0.01, 3), rep(0.8, 7)),
#'   significant = c(rep(TRUE, 3), rep(FALSE, 7))
#' )
#' plot_volcano(res)
plot_volcano <- function(
  data,
  lab_col = "tRNA",
  padj_cutoff = 0.05,
  max_overlaps = 20,
  point_size = 1.5,
  label_size = 3,
  sig_color = "#D55E00",
  nonsig_color = "grey60"
) {
  rlang::check_installed("ggrepel", reason = "to label significant points.")

  p <- ggplot(data, aes(x = log2FoldChange, y = -log10(pvalue))) +
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
    geom_hline(
      yintercept = -log10(padj_cutoff),
      linetype = "dashed",
      color = "grey40"
    ) +
    labs(
      x = "log2 Fold Change",
      y = "-log10(p-value)"
    ) +
    cowplot::theme_cowplot() +
    theme(legend.position = "none")

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
#'   `tRNA`, `log2FoldChange`, and `padj` columns.
#' @param charging_diffs A tibble from [compute_charging_diffs()] with
#'   at least `tRNA` and `diff` columns.
#' @param lab_col Column name (string) used for point labels. Default
#'   `"tRNA"`.
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
#'   tRNA = paste0("tRNA-", 1:6),
#'   log2FoldChange = c(1, -1, 0.5, -0.5, 2, -2),
#'   padj = c(0.01, 0.02, 0.5, 0.6, 0.001, 0.003)
#' )
#' charging_diffs <- tibble::tibble(
#'   tRNA = paste0("tRNA-", 1:6),
#'   diff = c(0.1, -0.1, 0.05, -0.05, -0.2, 0.15),
#'   se_diff = rep(0.03, 6)
#' )
#' plot_abundance_charging(deseq_res, charging_diffs)
plot_abundance_charging <- function(
  deseq_res,
  charging_diffs,
  lab_col = "tRNA",
  padj_cutoff = 0.05,
  max_overlaps = 20,
  point_size = 2,
  label_size = 3
) {
  rlang::check_installed("ggrepel", reason = "to label significant points.")

  data <- dplyr::inner_join(deseq_res, charging_diffs, by = "tRNA")

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
      x = "log2 Fold Change (abundance)",
      y = "Charging ratio difference",
      color = "Abundance / Charging"
    ) +
    cowplot::theme_cowplot() +
    theme(legend.position = "bottom")

  p
}

#' Plot per-tRNA charging ratio differences.
#'
#' Create a dot plot with error bars showing the difference in charging
#' ratio between two conditions for each tRNA. Consumes the tibble
#' returned by [compute_charging_diffs()].
#'
#' @param data A tibble from [compute_charging_diffs()] with at least
#'   `tRNA` (factor), `diff`, and `se_diff` columns.
#' @param point_size Numeric size for [ggplot2::geom_point()]. Default
#'   `2.5`.
#'
#' @return A ggplot object.
#'
#' @export
#'
#' @examples
#' df <- tibble::tibble(
#'   tRNA = forcats::fct_inorder(paste0("tRNA-", 1:5)),
#'   diff = c(-0.1, -0.05, 0.02, 0.08, 0.15),
#'   se_diff = rep(0.03, 5)
#' )
#' plot_charging_diffs(df)
plot_charging_diffs <- function(data, point_size = 2.5) {
  ggplot(data, aes(x = diff, y = tRNA)) +
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
