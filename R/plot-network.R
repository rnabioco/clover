# Network visualization functions -----------------------------------------------

#' Build a network from odds ratio data.
#'
#' Construct a tidygraph network from pairwise odds ratio or ROR data.
#' Each node is a position (or isodecoder) and each edge represents a
#' significant co-occurrence relationship. Requires the tidygraph and
#' igraph packages.
#'
#' @param data A tibble with columns `pos1`, `pos2`, and a value
#'   column specified by `value_col`. Optionally includes `p_adj`.
#' @param value_col Column name (string) for edge values. Default
#'   `"ror"`.
#' @param min_weight Minimum absolute value of `value_col` for an edge
#'   to be included. Default `0`.
#'
#' @return A `tbl_graph` object with node metrics `degree` and
#'   `betweenness`, or `NULL` if no edges remain after filtering.
#'
#' @export
#'
#' @examples
#' \donttest{
#' df <- tibble::tibble(
#'   pos1 = c("20", "34", "20"),
#'   pos2 = c("34", "58", "58"),
#'   ror = c(1.5, -0.8, 0.3)
#' )
#' build_or_network(df)
#' }
build_or_network <- function(data, value_col = "ror", min_weight = 0) {
  rlang::check_installed("tidygraph", reason = "to build tRNA networks.")
  rlang::check_installed("igraph", reason = "to compute network metrics.")

  edges <- data |>
    dplyr::filter(abs(.data[[value_col]]) > min_weight) |>
    dplyr::mutate(
      from = as.character(pmin(.data$pos1, .data$pos2)),
      to = as.character(pmax(.data$pos1, .data$pos2)),
      edge_type = dplyr::case_when(
        .data[[value_col]] > 0 ~ "co-occurring",
        .data[[value_col]] < 0 ~ "exclusive",
        .default = "neutral"
      ),
      weight = abs(.data[[value_col]])
    ) |>
    dplyr::select("from", "to", "edge_type", "weight")

  if (nrow(edges) == 0) {
    cli_warn("No edges remain after filtering.")
    return(NULL)
  }

  nodes <- tibble::tibble(
    position = unique(c(edges$from, edges$to))
  )

  graph <- tidygraph::tbl_graph(
    nodes = nodes,
    edges = edges,
    directed = FALSE
  )

  graph |>
    tidygraph::activate(nodes) |>
    dplyr::mutate(
      degree = tidygraph::centrality_degree(),
      betweenness = tidygraph::centrality_betweenness()
    )
}

#' Plot a network as an arc diagram.
#'
#' Create a circular arc diagram from a tidygraph network built by
#' [build_or_network()]. Requires the ggraph package.
#'
#' @param graph A `tbl_graph` object from [build_or_network()].
#' @param title Plot title. Default `"tRNA modification network"`.
#' @param co_color Color for co-occurring edges. Default
#'   `"steelblue"`.
#' @param ex_color Color for exclusive edges. Default `"tomato"`.
#'
#' @return A ggplot object, or `NULL` if `graph` is `NULL`.
#'
#' @export
#'
#' @examples
#' \donttest{
#' df <- tibble::tibble(
#'   pos1 = c("20", "34", "20"),
#'   pos2 = c("34", "58", "58"),
#'   ror = c(1.5, -0.8, 0.3)
#' )
#' graph <- build_or_network(df)
#' plot_arc_diagram(graph)
#' }
plot_arc_diagram <- function(
  graph,
  title = "tRNA modification network",
  co_color = "steelblue",
  ex_color = "tomato"
) {
  rlang::check_installed("ggraph", reason = "to draw arc diagrams.")

  if (is.null(graph)) {
    return(NULL)
  }

  ggraph::ggraph(graph, layout = "linear", circular = TRUE) +
    ggraph::geom_edge_arc(
      aes(
        edge_width = .data$weight,
        edge_alpha = .data$weight,
        edge_colour = .data$edge_type
      ),
      strength = 0.5
    ) +
    ggraph::scale_edge_colour_manual(
      values = c("co-occurring" = co_color, "exclusive" = ex_color),
      name = "Edge type"
    ) +
    ggraph::geom_node_point(aes(size = .data$degree), color = "grey30") +
    ggraph::geom_node_text(
      aes(label = .data$position),
      repel = TRUE,
      size = 3
    ) +
    ggraph::scale_edge_width(range = c(0.5, 3), name = "|ROR|") +
    ggraph::scale_edge_alpha(range = c(0.3, 0.9), guide = "none") +
    scale_size_continuous(range = c(2, 6), name = "Degree") +
    labs(title = title) +
    theme_void() +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
      legend.position = "bottom"
    )
}
