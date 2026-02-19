# tRNA secondary structure SVG visualization -----------------------------------

#' List organisms with bundled tRNA structure SVGs
#'
#' Returns the names of organisms for which tRNA cloverleaf structure
#' SVGs are bundled with the package. These organisms can be used
#' with [plot_tRNA_structure()].
#'
#' @return A character vector of organism names.
#'
#' @export
#'
#' @examples
#' structure_organisms()
structure_organisms <- function() {
  structures_dir <- system.file(
    "extdata",
    "structures",
    package = "clover"
  )
  if (structures_dir == "") {
    return(character(0))
  }
  dirs <- list.dirs(structures_dir, recursive = FALSE, full.names = FALSE)
  gsub("_", " ", dirs)
}

#' List available tRNA structures for an organism
#'
#' Returns the names of tRNAs for which cloverleaf structure SVGs
#' are bundled with the package for the given organism.
#'
#' @param organism Character string specifying the organism name
#'   (e.g., `"Escherichia coli"`). Use [structure_organisms()] to
#'   list available organisms.
#'
#' @return A character vector of tRNA names.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' structure_trnas("Escherichia coli")
#' }
structure_trnas <- function(organism) {
  org_dir <- structure_org_dir(organism)
  svg_files <- list.files(org_dir, pattern = "\\.svg$")
  tools::file_path_sans_ext(svg_files)
}

#' Plot tRNA secondary structure with modifications and linkages
#'
#' Reads a bundled tRNA cloverleaf SVG and overlays modification
#' highlights and circuit linkage arcs. Modifications are shown as
#' colored circles behind nucleotide letters; linkages are drawn as
#' Bezier curve arcs between position pairs.
#'
#' @param trna Character string identifying the tRNA
#'   (e.g., `"tRNA-Ala-GGC"`). Use [structure_trnas()] to list
#'   available tRNAs.
#' @param organism Character string specifying the organism name
#'   (e.g., `"Escherichia coli"`).
#' @param modifications A tibble with columns `pos` (1-based
#'   position in the tRNA sequence) and `mod1` (short modification
#'   name, e.g., `"m1A"`). Output of [modomics_mods()] works
#'   directly after filtering to the tRNA of interest.
#' @param linkages A tibble with columns `pos1`, `pos2`, and
#'   optionally `value` (e.g., log odds ratio) for coloring arcs.
#'   Output of [clean_odds_ratios()] works directly.
#' @param output Path for the output SVG file. If `NULL` (default),
#'   writes to a temporary file.
#' @param mod_palette Named character vector of colors keyed by
#'   modification short name. If `NULL`, uses a default palette.
#' @param linkage_palette Character vector of length 2 giving the
#'   low and high colors for the linkage value gradient. Default
#'   `c("#0072B2", "#D55E00")` (blue to vermillion).
#'
#' @return The path to the annotated SVG file (invisibly).
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Base structure only
#' plot_tRNA_structure("tRNA-Ala-GGC", "Escherichia coli")
#'
#' # With MODOMICS modifications
#' fa <- clover_example("ecoli/validated.fa.gz")
#' mods <- modomics_mods(fa, "Escherichia coli")
#' plot_tRNA_structure(
#'   "tRNA-Ala-GGC", "Escherichia coli",
#'   modifications = mods
#' )
#' }
plot_tRNA_structure <- function(
  trna,
  organism,
  modifications = NULL,
  linkages = NULL,
  output = NULL,
  mod_palette = NULL,
  linkage_palette = c("#0072B2", "#D55E00")
) {
  rlang::check_installed("jsonlite", reason = "to read structure metadata.")

  org_dir <- structure_org_dir(organism)

  svg_path <- file.path(org_dir, paste0(trna, ".svg"))
  json_path <- file.path(org_dir, paste0(trna, ".json"))

  if (!file.exists(svg_path)) {
    cli::cli_abort(c(
      "No structure SVG found for {.val {trna}}.",
      "i" = "Use {.fn structure_trnas} to list available tRNAs."
    ))
  }
  if (!file.exists(json_path)) {
    cli::cli_abort(
      "No position metadata found for {.val {trna}}."
    )
  }

  # Load base SVG and metadata
  svg_doc <- xml2::read_xml(svg_path)
  metadata <- jsonlite::fromJSON(json_path, simplifyVector = TRUE)
  nucs <- metadata$nucleotides

  if (is.null(output)) {
    output <- tempfile(fileext = ".svg")
  }

  # Add modification highlights
  if (!is.null(modifications)) {
    if (is.null(mod_palette)) {
      mod_palette <- default_mod_palette()
    }
    svg_doc <- add_mod_circles(svg_doc, nucs, modifications, mod_palette)
  }

  # Add linkage arcs
  if (!is.null(linkages)) {
    svg_doc <- add_linkage_arcs(svg_doc, nucs, linkages, linkage_palette)
  }

  # Add legend
  if (!is.null(modifications) || !is.null(linkages)) {
    svg_doc <- add_structure_legend(
      svg_doc,
      metadata,
      modifications,
      mod_palette,
      linkages,
      linkage_palette
    )
  }

  xml2::write_xml(svg_doc, output)
  cli::cli_inform("Wrote annotated SVG to {.path {output}}.")

  invisible(output)
}


# Internal helpers -------------------------------------------------------------

structure_org_dir <- function(organism) {
  org_fname <- gsub(" ", "_", organism)
  org_dir <- system.file(
    "extdata",
    "structures",
    org_fname,
    package = "clover"
  )
  if (org_dir == "") {
    cli::cli_abort(c(
      "No structure data found for {.val {organism}}.",
      "i" = "Use {.fn structure_organisms} to list available organisms."
    ))
  }
  org_dir
}

add_mod_circles <- function(svg_doc, nucs, modifications, palette) {
  svg_ns <- xml2::xml_ns(svg_doc)
  root <- xml2::xml_root(svg_doc)

  # Find the first child group to insert before (so circles are behind text)
  children <- xml2::xml_children(root)
  first_child <- if (length(children) > 0) children[[1]] else NULL

  # Create a group for modification circles
  mod_group <- xml2::xml_add_child(
    root,
    "g",
    id = "clover-modifications",
    .where = 0
  )

  for (i in seq_len(nrow(modifications))) {
    mod_pos <- modifications$pos[i]
    mod_name <- modifications$mod1[i]

    # Find matching nucleotide by position
    nuc_idx <- which(nucs$pos == mod_pos)
    if (length(nuc_idx) == 0) {
      next
    }

    nuc <- nucs[nuc_idx[1], ]
    color <- palette[mod_name]
    if (is.na(color) || is.null(color)) {
      color <- "#999999"
    }

    xml2::xml_add_child(
      mod_group,
      "circle",
      cx = as.character(nuc$x),
      cy = as.character(nuc$y),
      r = "6",
      fill = color,
      "fill-opacity" = "0.6",
      stroke = "none"
    )
  }

  svg_doc
}

add_linkage_arcs <- function(svg_doc, nucs, linkages, palette) {
  root <- xml2::xml_root(svg_doc)

  # Create a group for linkage arcs (behind modifications, behind text)
  arc_group <- xml2::xml_add_child(
    root,
    "g",
    id = "clover-linkages",
    .where = 0
  )

  has_value <- "value" %in% names(linkages)

  if (has_value) {
    vals <- linkages$value
    val_range <- range(vals, na.rm = TRUE)
  }

  for (i in seq_len(nrow(linkages))) {
    p1 <- linkages$pos1[i]
    p2 <- linkages$pos2[i]

    idx1 <- which(nucs$pos == p1)
    idx2 <- which(nucs$pos == p2)
    if (length(idx1) == 0 || length(idx2) == 0) {
      next
    }

    n1 <- nucs[idx1[1], ]
    n2 <- nucs[idx2[1], ]

    # Compute Bezier control point (offset perpendicular to midpoint)
    mx <- (n1$x + n2$x) / 2
    my <- (n1$y + n2$y) / 2
    dx <- n2$x - n1$x
    dy <- n2$y - n1$y
    dist <- sqrt(dx^2 + dy^2)

    if (dist < 1) {
      next
    }

    # Perpendicular offset (30% of distance)
    offset <- dist * 0.3
    cx <- mx - (dy / dist) * offset
    cy <- my + (dx / dist) * offset

    # Color based on value
    if (has_value && !is.na(linkages$value[i])) {
      color <- interpolate_color(
        linkages$value[i],
        val_range,
        palette
      )
    } else {
      color <- palette[1]
    }

    path_d <- sprintf(
      "M %.1f,%.1f Q %.1f,%.1f %.1f,%.1f",
      n1$x,
      n1$y,
      cx,
      cy,
      n2$x,
      n2$y
    )

    xml2::xml_add_child(
      arc_group,
      "path",
      d = path_d,
      fill = "none",
      stroke = color,
      "stroke-width" = "1.5",
      "stroke-opacity" = "0.7"
    )
  }

  svg_doc
}

add_structure_legend <- function(
  svg_doc,
  metadata,
  modifications,
  mod_palette,
  linkages,
  linkage_palette
) {
  root <- xml2::xml_root(svg_doc)
  svg_width <- metadata$width

  # Position legend to the right of the structure
  legend_x <- svg_width + 10
  legend_y <- 20

  legend_group <- xml2::xml_add_child(
    root,
    "g",
    id = "clover-legend",
    transform = sprintf("translate(%.0f, %.0f)", legend_x, legend_y)
  )

  y_offset <- 0

  # Modification legend
  if (!is.null(modifications) && !is.null(mod_palette)) {
    mod_types <- unique(modifications$mod1)
    mod_types <- mod_types[mod_types %in% names(mod_palette)]

    if (length(mod_types) > 0) {
      xml2::xml_add_child(
        legend_group,
        "text",
        x = "0",
        y = as.character(y_offset),
        "font-size" = "10",
        "font-weight" = "bold",
        "Modifications"
      )
      y_offset <- y_offset + 15

      for (mod in mod_types) {
        xml2::xml_add_child(
          legend_group,
          "circle",
          cx = "6",
          cy = as.character(y_offset - 3),
          r = "5",
          fill = mod_palette[mod],
          "fill-opacity" = "0.6"
        )
        xml2::xml_add_child(
          legend_group,
          "text",
          x = "16",
          y = as.character(y_offset),
          "font-size" = "9",
          mod
        )
        y_offset <- y_offset + 14
      }
    }
  }

  # Linkage legend
  if (!is.null(linkages)) {
    y_offset <- y_offset + 5
    xml2::xml_add_child(
      legend_group,
      "text",
      x = "0",
      y = as.character(y_offset),
      "font-size" = "10",
      "font-weight" = "bold",
      "Linkages"
    )
    y_offset <- y_offset + 15

    xml2::xml_add_child(
      legend_group,
      "line",
      x1 = "0",
      y1 = as.character(y_offset - 3),
      x2 = "20",
      y2 = as.character(y_offset - 3),
      stroke = linkage_palette[1],
      "stroke-width" = "1.5"
    )
    xml2::xml_add_child(
      legend_group,
      "text",
      x = "25",
      y = as.character(y_offset),
      "font-size" = "9",
      "co-occurrence"
    )
  }

  # Expand viewBox to accommodate legend
  viewbox <- xml2::xml_attr(root, "viewBox")
  if (!is.na(viewbox)) {
    parts <- as.numeric(strsplit(viewbox, "\\s+")[[1]])
    if (length(parts) == 4) {
      new_width <- parts[3] + 120
      xml2::xml_set_attr(
        root,
        "viewBox",
        paste(parts[1], parts[2], new_width, parts[4])
      )
    }
  }

  # Also update width attribute if present
  w_attr <- xml2::xml_attr(root, "width")
  if (!is.na(w_attr)) {
    w_val <- as.numeric(gsub("[^0-9.]", "", w_attr))
    if (!is.na(w_val)) {
      xml2::xml_set_attr(root, "width", as.character(w_val + 120))
    }
  }

  svg_doc
}

interpolate_color <- function(value, range, palette) {
  if (range[1] == range[2]) {
    return(palette[1])
  }
  t <- (value - range[1]) / (range[2] - range[1])
  t <- max(0, min(1, t))

  col1 <- grDevices::col2rgb(palette[1]) / 255
  col2 <- grDevices::col2rgb(palette[2]) / 255
  mixed <- col1 * (1 - t) + col2 * t

  grDevices::rgb(mixed[1], mixed[2], mixed[3])
}

default_mod_palette <- function() {
  c(
    "m1A" = "#E41A1C",
    "m5C" = "#377EB8",
    "m7G" = "#4DAF4A",
    "m1G" = "#984EA3",
    "m5U" = "#FF7F00",
    "D" = "#FFFF33",
    "pseudoU" = "#A65628",
    "I" = "#F781BF",
    "m2G" = "#66C2A5",
    "t6A" = "#FC8D62",
    "m22G" = "#8DA0CB",
    "Cm" = "#E78AC3",
    "Gm" = "#A6D854",
    "Um" = "#FFD92F",
    "Am" = "#E5C494"
  )
}
