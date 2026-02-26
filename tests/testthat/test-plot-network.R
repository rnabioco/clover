test_that("build_or_network returns tbl_graph", {
  skip_if_not_installed("tidygraph")
  skip_if_not_installed("igraph")

  df <- tibble::tibble(
    pos1 = c("20", "34", "20"),
    pos2 = c("34", "58", "58"),
    ror = c(1.5, -0.8, 0.3)
  )
  graph <- build_or_network(df)
  expect_s3_class(graph, "tbl_graph")
})

test_that("build_or_network filters by min_weight", {
  skip_if_not_installed("tidygraph")
  skip_if_not_installed("igraph")

  df <- tibble::tibble(
    pos1 = c("20", "34"),
    pos2 = c("34", "58"),
    ror = c(1.5, 0.1)
  )
  graph <- build_or_network(df, min_weight = 1.0)
  edges <- tidygraph::activate(graph, edges) |>
    tibble::as_tibble()
  expect_equal(nrow(edges), 1)
})

test_that("build_or_network returns NULL when no edges", {
  skip_if_not_installed("tidygraph")
  skip_if_not_installed("igraph")

  df <- tibble::tibble(
    pos1 = "20",
    pos2 = "34",
    ror = 0.1
  )
  expect_snapshot(graph <- build_or_network(df, min_weight = 5.0))
  expect_null(graph)
})

test_that("build_or_network classifies edge types", {
  skip_if_not_installed("tidygraph")
  skip_if_not_installed("igraph")

  df <- tibble::tibble(
    pos1 = c("20", "34"),
    pos2 = c("34", "58"),
    ror = c(1.5, -0.8)
  )
  graph <- build_or_network(df)
  edges <- tidygraph::activate(graph, edges) |>
    tibble::as_tibble()
  expect_equal(sort(edges$edge_type), c("co-occurring", "exclusive"))
})

test_that("plot_arc_diagram returns ggplot", {
  skip_if_not_installed("tidygraph")
  skip_if_not_installed("igraph")
  skip_if_not_installed("ggraph")

  df <- tibble::tibble(
    pos1 = c("20", "34"),
    pos2 = c("34", "58"),
    ror = c(1.5, -0.8)
  )
  graph <- build_or_network(df)
  p <- plot_arc_diagram(graph)
  expect_s3_class(p, "ggplot")
})

test_that("plot_arc_diagram returns NULL for NULL graph", {
  skip_if_not_installed("ggraph")

  expect_null(plot_arc_diagram(NULL))
})
