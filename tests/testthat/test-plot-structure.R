test_that("structure_organisms returns character vector", {
  result <- structure_organisms()
  expect_type(result, "character")
})

test_that("structure_trnas errors for unknown organism", {
  expect_snapshot(
    structure_trnas("Nonexistent organism"),
    error = TRUE
  )
})

test_that("plot_tRNA_structure errors for unknown organism", {
  expect_snapshot(
    plot_tRNA_structure("tRNA-Ala-GGC", "Nonexistent organism"),
    error = TRUE
  )
})

test_that("plot_tRNA_structure errors for unknown tRNA", {
  skip_if(
    length(structure_organisms()) == 0,
    "No bundled structure SVGs"
  )
  org <- structure_organisms()[1]
  expect_snapshot(
    plot_tRNA_structure("tRNA-Fake-XXX", org),
    error = TRUE
  )
})

test_that("default_mod_palette returns named character vector", {
  pal <- default_mod_palette()
  expect_type(pal, "character")
  expect_true(length(pal) > 0)
  expect_true(all(nchar(names(pal)) > 0))
})

test_that("interpolate_color returns valid color", {
  col <- interpolate_color(0.5, c(0, 1), c("#0072B2", "#D55E00"))
  expect_type(col, "character")
  expect_match(col, "^#[0-9A-Fa-f]{6}$")
})

test_that("interpolate_color handles equal range", {
  col <- interpolate_color(5, c(5, 5), c("#0072B2", "#D55E00"))
  expect_type(col, "character")
})

test_that("add_mod_circles annotates SVG with circles", {
  svg_text <- paste0(
    '<svg xmlns="http://www.w3.org/2000/svg" width="100" height="100">',
    '<text x="10" y="10">A</text>',
    '</svg>'
  )
  doc <- xml2::read_xml(svg_text)
  nucs <- data.frame(pos = 1L, base = "A", x = 10, y = 10)
  mods <- dplyr::tibble(pos = 1L, mod1 = "m1A")
  palette <- c(m1A = "#E41A1C")

  result <- add_mod_circles(doc, nucs, mods, palette)
  circles <- xml2::xml_find_all(result, ".//*[local-name()='circle']")
  expect_length(circles, 1)
})

test_that("add_linkage_arcs annotates SVG with paths", {
  svg_text <- paste0(
    '<svg xmlns="http://www.w3.org/2000/svg" width="200" height="200">',
    '<text x="10" y="10">A</text>',
    '<text x="50" y="50">G</text>',
    '</svg>'
  )
  doc <- xml2::read_xml(svg_text)
  nucs <- data.frame(
    pos = c(1L, 2L),
    base = c("A", "G"),
    x = c(10, 50),
    y = c(10, 50)
  )
  linkages <- dplyr::tibble(pos1 = 1L, pos2 = 2L, value = 1.5)
  palette <- c("#0072B2", "#D55E00")

  result <- add_linkage_arcs(doc, nucs, linkages, palette)
  paths <- xml2::xml_find_all(result, ".//*[local-name()='path']")
  expect_length(paths, 1)
})

test_that("add_linkage_arcs skips missing positions", {
  svg_text <- paste0(
    '<svg xmlns="http://www.w3.org/2000/svg" width="100" height="100">',
    '<text x="10" y="10">A</text>',
    '</svg>'
  )
  doc <- xml2::read_xml(svg_text)
  nucs <- data.frame(pos = 1L, base = "A", x = 10, y = 10)
  linkages <- dplyr::tibble(pos1 = 1L, pos2 = 99L)
  palette <- c("#0072B2", "#D55E00")

  result <- add_linkage_arcs(doc, nucs, linkages, palette)
  paths <- xml2::xml_find_all(result, ".//*[local-name()='path']")
  expect_length(paths, 0)
})
