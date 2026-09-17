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

test_that("plot_tRNA_structure auto-renames log_odds_ratio to value", {
  skip_if(
    length(structure_organisms()) == 0,
    "No bundled structure data"
  )
  org <- structure_organisms()[1]
  trna <- structure_trnas(org)[1]
  linkages <- dplyr::tibble(
    pos1 = 1L,
    pos2 = 10L,
    log_odds_ratio = 1.5
  )
  svg <- plot_tRNA_structure(trna, org, linkages = linkages)
  expect_true(file.exists(svg))
})

test_that("assign_arc_lanes assigns different lanes for overlapping arcs", {
  # Two arcs with similar angular spans should get different lanes
  arcs_info <- data.frame(
    arc_idx = 1:2,
    angle1 = c(0, 0.1),
    angle2 = c(1, 0.9)
  )
  lanes <- assign_arc_lanes(arcs_info)
  expect_length(lanes, 2)
  expect_false(lanes[1] == lanes[2])
})

test_that("assign_arc_lanes assigns same lane for non-overlapping arcs", {
  arcs_info <- data.frame(
    arc_idx = 1:2,
    angle1 = c(0, 2),
    angle2 = c(0.5, 2.5)
  )
  lanes <- assign_arc_lanes(arcs_info)
  expect_equal(lanes[1], lanes[2])
})

test_that("assign_arc_lanes handles empty input", {
  arcs_info <- data.frame(
    arc_idx = integer(),
    angle1 = numeric(),
    angle2 = numeric()
  )
  lanes <- assign_arc_lanes(arcs_info)
  expect_length(lanes, 0)
})

test_that("arc_stroke_width returns values in [1.0, 3.0]", {
  expect_equal(arc_stroke_width(0, c(0, 10)), 1.0)
  expect_equal(arc_stroke_width(10, c(0, 10)), 3.0)
  expect_equal(arc_stroke_width(5, c(0, 10)), 2.0)
})

test_that("arc_stroke_width handles equal range", {
  sw <- arc_stroke_width(5, c(5, 5))
  expect_equal(sw, 2.0)
})

test_that("add_linkage_arcs with bidirectional values produces paths", {
  svg_text <- paste0(
    '<svg xmlns="http://www.w3.org/2000/svg" width="200" height="200">',
    '<text x="10" y="10">A</text>',
    '<text x="50" y="50">G</text>',
    '<text x="90" y="10">C</text>',
    '</svg>'
  )
  doc <- xml2::read_xml(svg_text)
  nucs <- data.frame(
    pos = c(1L, 2L, 3L),
    base = c("A", "G", "C"),
    x = c(10, 50, 90),
    y = c(10, 50, 10)
  )
  linkages <- dplyr::tibble(
    pos1 = c(1L, 1L),
    pos2 = c(2L, 3L),
    value = c(-1.5, 2.0)
  )
  palette <- c("#0072B2", "#D55E00")

  result <- add_linkage_arcs(doc, nucs, linkages, palette)
  paths <- xml2::xml_find_all(result, ".//*[local-name()='path']")
  expect_length(paths, 2)

  # First arc (negative value) should use palette[1]
  expect_equal(xml2::xml_attr(paths[[1]], "stroke"), "#0072B2")
  # Second arc (positive value) should use palette[2]
  expect_equal(xml2::xml_attr(paths[[2]], "stroke"), "#D55E00")
})

test_that("structure_to_png converts SVG to PNG", {
  skip_if_not_installed("rsvg")
  skip_if(
    length(structure_organisms()) == 0,
    "No bundled structure SVGs"
  )

  org <- structure_organisms()[1]
  trna <- structure_trnas(org)[1]
  svg <- plot_tRNA_structure(trna, org)

  png <- structure_to_png(svg)
  expect_true(file.exists(png))
  expect_match(png, "\\.png$")
  expect_gt(file.info(png)$size, 0)
})

test_that("structure_to_png errors for missing file", {
  skip_if_not_installed("rsvg")
  expect_snapshot(
    structure_to_png("nonexistent.svg"),
    error = TRUE
  )
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

test_that("add_end_labels adds amino acid label and line", {
  svg_text <- paste0(
    '<svg xmlns="http://www.w3.org/2000/svg" width="100" height="100">',
    '<text x="10" y="10">G</text>',
    '</svg>'
  )
  doc <- xml2::read_xml(svg_text)
  nucs <- data.frame(pos = 1L, base = "G", x = 10, y = 10)
  metadata <- list(trna_name = "tRNA-Glu-TTC")

  result <- add_end_labels(doc, nucs, metadata)
  group <- xml2::xml_find_first(result, ".//*[@id='clover-end-labels']")
  expect_false(is.na(group))

  lines <- xml2::xml_find_all(group, ".//*[local-name()='line']")
  expect_length(lines, 1)

  texts <- xml2::xml_find_all(group, ".//*[local-name()='text']")
  expect_length(texts, 1)
  expect_equal(xml2::xml_text(texts[[1]]), "Glu")
})

test_that("add_position_markers adds markers every 10 nt", {
  svg_text <- paste0(
    '<svg xmlns="http://www.w3.org/2000/svg" width="200" height="200">',
    '<text x="10" y="10">A</text>',
    '</svg>'
  )
  doc <- xml2::read_xml(svg_text)
  nucs <- data.frame(
    pos = 1:25,
    base = rep("A", 25),
    x = seq(10, 250, by = 10),
    y = rep(50, 25)
  )

  result <- add_position_markers(doc, nucs)
  group <- xml2::xml_find_first(result, ".//*[@id='clover-position-markers']")
  expect_false(is.na(group))

  texts <- xml2::xml_find_all(group, ".//*[local-name()='text']")
  expect_length(texts, 2)
  expect_equal(xml2::xml_text(texts[[1]]), "10")
  expect_equal(xml2::xml_text(texts[[2]]), "20")
})

test_that("structure_html returns html class", {
  skip_if_not_installed("htmltools")
  svg_file <- withr::local_tempfile(fileext = ".svg")
  writeLines(
    '<svg xmlns="http://www.w3.org/2000/svg"><text>A</text></svg>',
    svg_file
  )

  result <- structure_html(svg_file)
  expect_s3_class(result, "html")
  expect_match(as.character(result), "text-align: center")
})

test_that("structure_html errors on missing file", {
  skip_if_not_installed("htmltools")
  expect_snapshot(
    structure_html("nonexistent.svg"),
    error = TRUE
  )
})

# Sprinzl coordinate conversion ------------------------------------------------

test_that("convert_sprinzl_positions maps labels to positions", {
  trna_coords <- dplyr::tibble(
    sprinzl_label = c("1", "2", "34", "35", "36"),
    pos = c(1L, 2L, 30L, 31L, 32L)
  )
  df <- dplyr::tibble(pos = c("34", "35", "36"), mod1 = c("m1A", "m5C", "D"))
  result <- convert_sprinzl_positions(df, "pos", trna_coords)
  expect_equal(result$pos, c(30L, 31L, 32L))
  expect_equal(result$mod1, c("m1A", "m5C", "D"))
})

test_that("convert_sprinzl_positions coerces numeric input", {
  trna_coords <- dplyr::tibble(
    sprinzl_label = c("34", "35"),
    pos = c(30L, 31L)
  )
  df <- dplyr::tibble(pos = c(34, 35), mod1 = c("m1A", "m5C"))
  result <- convert_sprinzl_positions(df, "pos", trna_coords)
  expect_equal(result$pos, c(30L, 31L))
})

test_that("convert_sprinzl_positions warns on unmatched and drops rows", {
  trna_coords <- dplyr::tibble(
    sprinzl_label = c("1", "2"),
    pos = c(1L, 2L)
  )
  df <- dplyr::tibble(pos = c("1", "99"), mod1 = c("m1A", "m5C"))
  expect_snapshot(
    result <- convert_sprinzl_positions(df, "pos", trna_coords)
  )
  expect_equal(nrow(result), 1)
  expect_equal(result$pos, 1L)
})

test_that("convert_sprinzl_positions converts two columns for linkages", {
  trna_coords <- dplyr::tibble(
    sprinzl_label = c("34", "35", "36"),
    pos = c(30L, 31L, 32L)
  )
  df <- dplyr::tibble(
    pos1 = c("34", "35"),
    pos2 = c("36", "34"),
    value = c(1.5, -0.5)
  )
  result <- convert_sprinzl_positions(df, c("pos1", "pos2"), trna_coords)
  expect_equal(result$pos1, c(30L, 31L))
  expect_equal(result$pos2, c(32L, 30L))
})

test_that("plot_tRNA_structure errors when tRNA not in sprinzl_coords", {
  skip_if(
    length(structure_organisms()) == 0,
    "No bundled structure SVGs"
  )

  org <- structure_organisms()[1]
  trna <- structure_trnas(org)[1]
  fake_coords <- dplyr::tibble(
    trna_id = "nuc-tRNA-Fake-AAA-1-1",
    pos = 1L,
    sprinzl_label = "1"
  )
  expect_snapshot(
    plot_tRNA_structure(trna, org, sprinzl_coords = fake_coords),
    error = TRUE
  )
})

test_that("plot_tRNA_structure converts sprinzl coords with real data", {
  skip_if(
    length(structure_organisms()) == 0,
    "No bundled structure SVGs"
  )
  coords_path <- system.file(
    "extdata",
    "sprinzl",
    "ecoliK12_global_coords.tsv.gz",
    package = "clover"
  )
  skip_if(coords_path == "", "No bundled sprinzl coords")

  coords <- read_sprinzl_coords(coords_path)
  org <- "Escherichia coli"
  trna <- "tRNA-Glu-TTC"

  mods <- dplyr::tibble(pos = c("34", "35"), mod1 = c("m1A", "m5C"))
  svg <- plot_tRNA_structure(
    trna,
    org,
    modifications = mods,
    sprinzl_coords = coords
  )
  expect_true(file.exists(svg))
})

test_that("plot_tRNA_structure respects position_markers = FALSE", {
  skip_if(
    length(structure_organisms()) == 0,
    "No bundled structure SVGs"
  )

  org <- structure_organisms()[1]
  trna <- structure_trnas(org)[1]

  svg_path <- plot_tRNA_structure(trna, org, position_markers = FALSE)
  doc <- xml2::read_xml(svg_path)
  markers <- xml2::xml_find_first(doc, ".//*[@id='clover-position-markers']")
  expect_true(is.na(markers))

  # End labels should still be present
  end_labels <- xml2::xml_find_first(doc, ".//*[@id='clover-end-labels']")
  expect_false(is.na(end_labels))
})
