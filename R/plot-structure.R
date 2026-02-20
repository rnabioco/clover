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
#' highlights, outline circles, and circuit linkage arcs.
#' Modifications are shown as colored filled circles behind
#' nucleotide letters; outlines are shown as colored circle borders;
#' linkages are drawn as Bezier curve arcs between position pairs.
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
#' @param outlines A tibble with columns `pos` (1-based position)
#'   and `group` (category name for palette lookup). Draws circle
#'   outlines (stroke only, no fill) around each nucleotide.
#' @param linkages A tibble with columns `pos1`, `pos2`, and
#'   optionally `value` (e.g., log odds ratio) for coloring arcs.
#'   Output of [clean_odds_ratios()] works directly.
#' @param output Path for the output SVG file. If `NULL` (default),
#'   writes to a temporary file.
#' @param mod_palette Named character vector of colors keyed by
#'   modification short name. If `NULL`, uses a default palette.
#' @param outline_palette Named character vector of colors keyed by
#'   outline group name. If `NULL`, uses `"#333333"` for all.
#' @param text_colors A tibble with columns `pos` (1-based position)
#'   and `color` (hex color string). Changes the nucleotide letter
#'   color at specified positions. Unspecified positions keep the
#'   default color.
#' @param position_markers Logical; if `TRUE` (default), draw
#'   small grey position numbers every 10 nucleotides around the
#'   cloverleaf to help orient readers.
#' @param linkage_palette Character vector of length 2 giving the
#'   colors for negative (exclusive) and positive (co-occurring)
#'   linkage values. Default `c("#0072B2", "#D55E00")` (blue for
#'   exclusive, vermillion for co-occurring). Stroke width encodes
#'   the magnitude of the value.
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
  outlines = NULL,
  linkages = NULL,
  output = NULL,
  mod_palette = NULL,
  outline_palette = NULL,
  text_colors = NULL,
  position_markers = TRUE,
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

  # Restyle base SVG: black nucleotide text, grey base-pair lines
  svg_doc <- restyle_base_svg(svg_doc)

  # Add 3' amino acid label
  svg_doc <- add_end_labels(svg_doc, nucs, metadata)

  # Add position markers every 10 nt
  if (position_markers) {
    svg_doc <- add_position_markers(svg_doc, nucs)
  }

  # Add modification highlights (filled circles behind text)
  if (!is.null(modifications)) {
    if (is.null(mod_palette)) {
      mod_palette <- default_mod_palette()
    }
    svg_doc <- add_mod_circles(svg_doc, nucs, modifications, mod_palette)
  }

  # Add outline circles (stroke-only circles on top of fills, behind text)
  if (!is.null(outlines)) {
    svg_doc <- add_outline_circles(
      svg_doc,
      nucs,
      outlines,
      outline_palette
    )
  }

  # Recolor nucleotide text at specified positions
  if (!is.null(text_colors)) {
    svg_doc <- recolor_text(svg_doc, nucs, text_colors)
  }

  # Add linkage arcs
  if (!is.null(linkages)) {
    svg_doc <- add_linkage_arcs(svg_doc, nucs, linkages, linkage_palette)
  }

  # Add legend
  if (!is.null(modifications) || !is.null(outlines) || !is.null(linkages)) {
    svg_doc <- add_structure_legend(
      svg_doc,
      metadata,
      modifications,
      mod_palette,
      outlines,
      outline_palette,
      linkages,
      linkage_palette
    )
  }

  xml2::write_xml(svg_doc, output)

  invisible(output)
}

#' Convert a tRNA structure SVG to PNG
#'
#' Renders an SVG file (typically produced by [plot_tRNA_structure()])
#' to a PNG bitmap. Requires the
#' \href{https://cran.r-project.org/package=rsvg}{rsvg} package.
#'
#' @param svg_path Path to an SVG file, typically the return value of
#'   [plot_tRNA_structure()].
#' @param output Path for the output PNG file. If `NULL` (default),
#'   replaces the `.svg` extension with `.png`.
#' @param width Width of the output PNG in pixels. If `NULL`
#'   (default), uses the intrinsic SVG width.
#' @param height Height of the output PNG in pixels. If `NULL`
#'   (default), uses the intrinsic SVG height.
#'
#' @return The path to the PNG file (invisibly).
#'
#' @export
#'
#' @examples
#' \dontrun{
#' svg <- plot_tRNA_structure("tRNA-Glu-TTC", "Escherichia coli")
#' png <- structure_to_png(svg)
#' }
structure_to_png <- function(
  svg_path,
  output = NULL,
  width = NULL,
  height = NULL
) {
  rlang::check_installed("rsvg", reason = "to convert SVG to PNG.")

  if (!file.exists(svg_path)) {
    cli::cli_abort("SVG file not found: {.path {svg_path}}.")
  }

  if (is.null(output)) {
    output <- sub("\\.svg$", ".png", svg_path, ignore.case = TRUE)
    if (output == svg_path) {
      output <- paste0(svg_path, ".png")
    }
  }

  rsvg::rsvg_png(svg_path, file = output, width = width, height = height)

  invisible(output)
}


# Internal helpers -------------------------------------------------------------

# R2R SVGs use font-size 7.5 Helvetica. The text x/y attributes give the
# left baseline of the character. These offsets shift to the visual center
# of the uppercase letter (approximately half character-width right, half
# cap-height up).
nuc_x_offset <- 2.5
nuc_y_offset <- -2.7

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

restyle_base_svg <- function(svg_doc) {
  root <- xml2::xml_root(svg_doc)
  ns <- xml2::xml_ns(svg_doc)

  # Change all nucleotide tspan text to black (R2R uses #d90000 red)
  tspans <- xml2::xml_find_all(root, ".//d1:tspan", ns)
  for (ts in tspans) {
    fill <- xml2::xml_attr(ts, "fill")
    if (!is.na(fill) && fill == "#d90000") {
      xml2::xml_set_attr(ts, "fill", "#000000")
    }
  }

  # Change base-pair line paths to grey (R2R uses black #000000 with
  # stroke-width 1.44)
  paths <- xml2::xml_find_all(root, ".//d1:path", ns)
  for (p in paths) {
    sw <- xml2::xml_attr(p, "stroke-width")
    stroke <- xml2::xml_attr(p, "stroke")
    if (!is.na(sw) && !is.na(stroke) && sw == "1.44" && stroke == "#000000") {
      xml2::xml_set_attr(p, "stroke", "#999999")
    }
  }

  svg_doc
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
      cx = as.character(nuc$x + nuc_x_offset),
      cy = as.character(nuc$y + nuc_y_offset),
      r = "4.8",
      fill = color,
      "fill-opacity" = "0.6",
      stroke = "none"
    )
  }

  svg_doc
}

recolor_text <- function(svg_doc, nucs, text_colors) {
  root <- xml2::xml_root(svg_doc)

  # Find all tspan elements (nucleotide letters)
  tspans <- xml2::xml_find_all(root, ".//d1:tspan", xml2::xml_ns(svg_doc))

  for (i in seq_len(nrow(text_colors))) {
    tc_pos <- text_colors$pos[i]
    tc_color <- text_colors$color[i]

    nuc_idx <- which(nucs$pos == tc_pos)
    if (length(nuc_idx) == 0) {
      next
    }

    nuc <- nucs[nuc_idx[1], ]

    # Match tspan by x coordinate (R2R sets x on both <text> and <tspan>)
    for (ts in tspans) {
      tx <- as.numeric(xml2::xml_attr(ts, "x"))
      ty <- as.numeric(xml2::xml_attr(ts, "y"))
      if (
        !is.na(tx) &&
          !is.na(ty) &&
          abs(tx - nuc$x) < 0.01 &&
          abs(ty - nuc$y) < 0.01
      ) {
        xml2::xml_set_attr(ts, "fill", tc_color)
        break
      }
    }
  }

  svg_doc
}

add_outline_circles <- function(svg_doc, nucs, outlines, palette) {
  root <- xml2::xml_root(svg_doc)

  # Insert after modifications group (index 1) if it exists, else at 0
  mod_group <- xml2::xml_find_first(root, ".//g[@id='clover-modifications']")
  where <- if (is.na(mod_group)) 0L else 1L

  outline_group <- xml2::xml_add_child(
    root,
    "g",
    id = "clover-outlines",
    .where = where
  )

  for (i in seq_len(nrow(outlines))) {
    out_pos <- outlines$pos[i]
    out_group <- outlines$group[i]

    nuc_idx <- which(nucs$pos == out_pos)
    if (length(nuc_idx) == 0) {
      next
    }

    nuc <- nucs[nuc_idx[1], ]
    color <- if (!is.null(palette)) palette[out_group] else NA
    if (is.na(color)) {
      color <- "#333333"
    }

    xml2::xml_add_child(
      outline_group,
      "circle",
      cx = as.character(nuc$x + nuc_x_offset),
      cy = as.character(nuc$y + nuc_y_offset),
      r = "4.8",
      fill = "none",
      stroke = color,
      "stroke-width" = "1.2"
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

  # Compute centroid of all nucleotide visual centers
  centroid_x <- mean(nucs$x + nuc_x_offset)
  centroid_y <- mean(nucs$y + nuc_y_offset)

  base_offset <- 20

  # Resolve coordinates for each arc and build info for lane assignment
  arcs <- list()
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

    # Use visual centers for all geometry
    n1_cx <- n1$x + nuc_x_offset
    n1_cy <- n1$y + nuc_y_offset
    n2_cx <- n2$x + nuc_x_offset
    n2_cy <- n2$y + nuc_y_offset

    if (sqrt((n2_cx - n1_cx)^2 + (n2_cy - n1_cy)^2) < 1) {
      next
    }

    # Angular span relative to centroid (for lane assignment)
    angle1 <- atan2(n1_cy - centroid_y, n1_cx - centroid_x)
    angle2 <- atan2(n2_cy - centroid_y, n2_cx - centroid_x)

    arcs[[length(arcs) + 1]] <- list(
      idx = i,
      n1_cx = n1_cx,
      n1_cy = n1_cy,
      n2_cx = n2_cx,
      n2_cy = n2_cy,
      angle1 = angle1,
      angle2 = angle2
    )
  }

  if (length(arcs) == 0) {
    return(svg_doc)
  }

  # Assign lanes to avoid overlap
  arcs_info <- data.frame(
    arc_idx = seq_along(arcs),
    angle1 = vapply(arcs, \(a) a$angle1, numeric(1)),
    angle2 = vapply(arcs, \(a) a$angle2, numeric(1))
  )
  lanes <- assign_arc_lanes(arcs_info)

  # Compute abs value range for stroke width mapping
  if (has_value) {
    abs_vals <- abs(linkages$value[
      !is.na(linkages$value) &
        is.finite(linkages$value)
    ])
    abs_range <- if (length(abs_vals) > 0) range(abs_vals) else c(0, 0)
  }

  for (j in seq_along(arcs)) {
    a <- arcs[[j]]
    i <- a$idx
    lane <- lanes[j]

    # Midpoint of the two visual centers
    mx <- (a$n1_cx + a$n2_cx) / 2
    my <- (a$n1_cy + a$n2_cy) / 2

    # Vector from centroid to midpoint
    vx <- mx - centroid_x
    vy <- my - centroid_y
    vmag <- sqrt(vx^2 + vy^2)

    # Offset distance scaled by lane
    offset <- base_offset * (1.0 + (lane - 1) * 0.7)

    if (vmag > 0.01) {
      # Extend outward from centroid
      cx <- mx + (vx / vmag) * offset
      cy <- my + (vy / vmag) * offset
    } else {
      # Fallback: perpendicular offset when midpoint is at centroid
      dx <- a$n2_cx - a$n1_cx
      dy <- a$n2_cy - a$n1_cy
      dist <- sqrt(dx^2 + dy^2)
      cx <- mx - (dy / dist) * offset
      cy <- my + (dx / dist) * offset
    }

    # Color by sign of value
    if (has_value && !is.na(linkages$value[i])) {
      color <- if (linkages$value[i] < 0) palette[1] else palette[2]
      abs_val <- if (is.finite(linkages$value[i])) {
        abs(linkages$value[i])
      } else {
        abs_range[2]
      }
      sw <- arc_stroke_width(abs_val, abs_range)
    } else {
      color <- palette[2]
      sw <- 1.5
    }

    path_d <- sprintf(
      "M %.1f,%.1f Q %.1f,%.1f %.1f,%.1f",
      a$n1_cx,
      a$n1_cy,
      cx,
      cy,
      a$n2_cx,
      a$n2_cy
    )

    xml2::xml_add_child(
      arc_group,
      "path",
      d = path_d,
      fill = "none",
      stroke = color,
      "stroke-width" = as.character(round(sw, 1)),
      "stroke-opacity" = "0.7"
    )
  }

  svg_doc
}

assign_arc_lanes <- function(arcs_info) {
  n <- nrow(arcs_info)
  if (n == 0) {
    return(integer(0))
  }

  # Normalize angular spans to [start, end] where start < end on the circle
  spans <- lapply(seq_len(n), function(i) {
    a1 <- arcs_info$angle1[i]
    a2 <- arcs_info$angle2[i]
    # Ensure consistent ordering: smaller angular span
    diff <- (a2 - a1) %% (2 * pi)
    if (diff > pi) {
      list(start = a2, end = a1 + 2 * pi)
    } else {
      list(start = a1, end = a2)
    }
  })

  # Sort by span size (smallest first)
  span_sizes <- vapply(
    spans,
    function(s) (s$end - s$start) %% (2 * pi),
    numeric(1)
  )
  order_idx <- order(span_sizes)

  lanes <- integer(n)
  lane_assignments <- list() # list of lists, one per lane

  for (idx in order_idx) {
    assigned <- FALSE
    for (lane_num in seq_along(lane_assignments)) {
      conflict <- FALSE
      for (other_idx in lane_assignments[[lane_num]]) {
        if (angular_spans_overlap(spans[[idx]], spans[[other_idx]])) {
          conflict <- TRUE
          break
        }
      }
      if (!conflict) {
        lanes[idx] <- lane_num
        lane_assignments[[lane_num]] <- c(lane_assignments[[lane_num]], idx)
        assigned <- TRUE
        break
      }
    }
    if (!assigned) {
      new_lane <- length(lane_assignments) + 1
      lanes[idx] <- new_lane
      lane_assignments[[new_lane]] <- idx
    }
  }

  lanes
}

angular_spans_overlap <- function(arc_a, arc_b) {
  # Normalize both spans to start in [0, 2*pi)
  twopi <- 2 * pi
  a_start <- arc_a$start %% twopi
  a_end <- a_start + ((arc_a$end - arc_a$start) %% twopi)
  b_start <- arc_b$start %% twopi
  b_end <- b_start + ((arc_b$end - arc_b$start) %% twopi)

  # Check overlap on the line; also check with b shifted by 2*pi
  overlaps <- function(s1, e1, s2, e2) {
    s1 < e2 && s2 < e1
  }

  overlaps(a_start, a_end, b_start, b_end) ||
    overlaps(a_start, a_end, b_start + twopi, b_end + twopi) ||
    overlaps(a_start + twopi, a_end + twopi, b_start, b_end)
}

arc_stroke_width <- function(abs_value, abs_range) {
  min_width <- 1.0
  max_width <- 3.0
  if (abs_range[1] == abs_range[2]) {
    return((min_width + max_width) / 2)
  }
  t <- (abs_value - abs_range[1]) / (abs_range[2] - abs_range[1])
  t <- max(0, min(1, t))
  min_width + t * (max_width - min_width)
}

add_structure_legend <- function(
  svg_doc,
  metadata,
  modifications,
  mod_palette,
  outlines,
  outline_palette,
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
        color <- mod_palette[mod]
        if (is.na(color) || is.null(color)) {
          color <- "#999999"
        }
        xml2::xml_add_child(
          legend_group,
          "circle",
          cx = "6",
          cy = as.character(y_offset - 3),
          r = "5",
          fill = color,
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

  # Outline legend
  if (!is.null(outlines) && !is.null(outline_palette)) {
    out_types <- unique(outlines$group)
    out_types <- out_types[out_types %in% names(outline_palette)]

    if (length(out_types) > 0) {
      y_offset <- y_offset + 5
      xml2::xml_add_child(
        legend_group,
        "text",
        x = "0",
        y = as.character(y_offset),
        "font-size" = "10",
        "font-weight" = "bold",
        "Outlines"
      )
      y_offset <- y_offset + 15

      for (out in out_types) {
        xml2::xml_add_child(
          legend_group,
          "circle",
          cx = "6",
          cy = as.character(y_offset - 3),
          r = "5",
          fill = "none",
          stroke = outline_palette[out],
          "stroke-width" = "1.2"
        )
        xml2::xml_add_child(
          legend_group,
          "text",
          x = "16",
          y = as.character(y_offset),
          "font-size" = "9",
          out
        )
        y_offset <- y_offset + 14
      }
    }
  }

  # Linkage legend
  if (!is.null(linkages)) {
    has_value <- "value" %in% names(linkages)
    has_neg <- has_value && any(linkages$value < 0, na.rm = TRUE)
    has_pos <- has_value && any(linkages$value >= 0, na.rm = TRUE)

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

    if (has_value && has_neg && has_pos) {
      # Bidirectional: show both entries
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
        "Exclusive"
      )
      y_offset <- y_offset + 14

      xml2::xml_add_child(
        legend_group,
        "line",
        x1 = "0",
        y1 = as.character(y_offset - 3),
        x2 = "20",
        y2 = as.character(y_offset - 3),
        stroke = linkage_palette[2],
        "stroke-width" = "1.5"
      )
      xml2::xml_add_child(
        legend_group,
        "text",
        x = "25",
        y = as.character(y_offset),
        "font-size" = "9",
        "Co-occurring"
      )
    } else if (has_value && has_neg) {
      # All negative
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
        "Exclusive"
      )
    } else if (has_value && has_pos) {
      # All positive
      xml2::xml_add_child(
        legend_group,
        "line",
        x1 = "0",
        y1 = as.character(y_offset - 3),
        x2 = "20",
        y2 = as.character(y_offset - 3),
        stroke = linkage_palette[2],
        "stroke-width" = "1.5"
      )
      xml2::xml_add_child(
        legend_group,
        "text",
        x = "25",
        y = as.character(y_offset),
        "font-size" = "9",
        "Co-occurring"
      )
    } else {
      # No value column
      xml2::xml_add_child(
        legend_group,
        "line",
        x1 = "0",
        y1 = as.character(y_offset - 3),
        x2 = "20",
        y2 = as.character(y_offset - 3),
        stroke = linkage_palette[2],
        "stroke-width" = "1.5"
      )
      xml2::xml_add_child(
        legend_group,
        "text",
        x = "25",
        y = as.character(y_offset),
        "font-size" = "9",
        "Linkage"
      )
    }
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

add_end_labels <- function(svg_doc, nucs, metadata) {
  root <- xml2::xml_root(svg_doc)

  # Extract amino acid from tRNA name (e.g., "tRNA-Glu-TTC" -> "Glu")
  trna_name <- metadata$trna_name
  aa <- strsplit(trna_name, "-")[[1]][2]

  # Remove R2R's amino acid indicator text (single letter beyond last

  # nucleotide that doesn't match any position in the metadata)
  remove_extra_aa_text(svg_doc, nucs)

  # Get last nucleotide (3' end)
  last_nuc <- nucs[nrow(nucs), ]
  last_cx <- last_nuc$x + nuc_x_offset
  last_cy <- last_nuc$y + nuc_y_offset

  # Compute centroid of all nucleotide visual centers
  centroid_x <- mean(nucs$x + nuc_x_offset)
  centroid_y <- mean(nucs$y + nuc_y_offset)

  # Direction from centroid toward the last nucleotide
  dx <- last_cx - centroid_x
  dy <- last_cy - centroid_y
  mag <- sqrt(dx^2 + dy^2)

  if (mag > 0.01) {
    dir_x <- dx / mag
    dir_y <- dy / mag
  } else {
    dir_x <- 0
    dir_y <- -1
  }

  # Label position: offset outward from the nucleotide
  offset <- 15
  label_x <- last_cx + dir_x * offset
  label_y <- last_cy + dir_y * offset

  label_group <- xml2::xml_add_child(
    root,
    "g",
    id = "clover-end-labels"
  )

  # Line from nucleotide center to label
  xml2::xml_add_child(
    label_group,
    "line",
    x1 = as.character(last_cx),
    y1 = as.character(last_cy),
    x2 = as.character(label_x),
    y2 = as.character(label_y),
    stroke = "black",
    "stroke-width" = "0.5"
  )

  # Rounded rectangle behind the label
  # Approximate text width based on character count
  pad_x <- 3
  pad_y <- 2
  text_width <- nchar(aa) * 4.5
  rect_w <- text_width + 2 * pad_x
  rect_h <- 7.5 + 2 * pad_y

  xml2::xml_add_child(
    label_group,
    "rect",
    x = as.character(label_x - rect_w / 2),
    y = as.character(label_y - rect_h / 2),
    width = as.character(rect_w),
    height = as.character(rect_h),
    rx = "3",
    ry = "3",
    fill = "white",
    stroke = "black",
    "stroke-width" = "0.8"
  )

  # Amino acid text label
  xml2::xml_add_child(
    label_group,
    "text",
    x = as.character(label_x),
    y = as.character(label_y),
    "font-size" = "7.5",
    "font-family" = "Helvetica, Arial, sans-serif",
    "text-anchor" = "middle",
    "dominant-baseline" = "central",
    fill = "black",
    aa
  )

  svg_doc
}

remove_extra_aa_text <- function(svg_doc, nucs) {
  root <- xml2::xml_root(svg_doc)
  tspans <- xml2::xml_find_all(root, ".//d1:tspan", xml2::xml_ns(svg_doc))

  # Build lookup of all known nucleotide coordinates
  nuc_coords <- data.frame(
    x = nucs$x,
    y = nucs$y
  )

  for (ts in tspans) {
    tx <- as.numeric(xml2::xml_attr(ts, "x"))
    ty <- as.numeric(xml2::xml_attr(ts, "y"))
    text_content <- xml2::xml_text(ts)

    if (is.na(tx) || is.na(ty)) {
      next
    }
    # Only look for single-character nucleotide letters
    if (nchar(text_content) != 1) {
      next
    }
    if (!text_content %in% c("A", "C", "G", "U", "T")) {
      next
    }

    # Check if this tspan matches any known nucleotide position
    dists <- sqrt((nuc_coords$x - tx)^2 + (nuc_coords$y - ty)^2)
    if (min(dists) > 0.1) {
      # This text element doesn't match any nucleotide — it's the extra
      # amino acid indicator from R2R. Hide it.
      parent <- xml2::xml_parent(ts)
      xml2::xml_set_attr(parent, "display", "none")
    }
  }
}

add_position_markers <- function(svg_doc, nucs) {
  root <- xml2::xml_root(svg_doc)

  marker_group <- xml2::xml_add_child(
    root,
    "g",
    id = "clover-position-markers"
  )

  max_pos <- max(nucs$pos)
  marker_positions <- seq(10, max_pos, by = 10)

  offset <- 10

  # Precompute visual centers for all nucleotides
  all_cx <- nucs$x + nuc_x_offset
  all_cy <- nucs$y + nuc_y_offset

  # Try 12 evenly-spaced candidate directions and pick the one whose

  # endpoint is furthest from any other nucleotide
  n_candidates <- 12
  angles <- seq(0, 2 * pi, length.out = n_candidates + 1)[-(n_candidates + 1)]

  for (pos in marker_positions) {
    nuc_idx <- which(nucs$pos == pos)
    if (length(nuc_idx) == 0) {
      next
    }

    nuc <- nucs[nuc_idx[1], ]
    nuc_cx <- nuc$x + nuc_x_offset
    nuc_cy <- nuc$y + nuc_y_offset

    best_min_dist <- -1
    best_x <- nuc_cx
    best_y <- nuc_cy - offset

    for (angle in angles) {
      cand_x <- nuc_cx + cos(angle) * offset
      cand_y <- nuc_cy + sin(angle) * offset

      dists <- sqrt((all_cx - cand_x)^2 + (all_cy - cand_y)^2)
      min_dist <- min(dists)

      if (min_dist > best_min_dist) {
        best_min_dist <- min_dist
        best_x <- cand_x
        best_y <- cand_y
      }
    }

    xml2::xml_add_child(
      marker_group,
      "text",
      x = as.character(best_x),
      y = as.character(best_y),
      "font-size" = "5.5",
      "font-family" = "Helvetica, Arial, sans-serif",
      "text-anchor" = "middle",
      "dominant-baseline" = "central",
      fill = "#666666",
      as.character(pos)
    )
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

#' Default modification color palette
#'
#' Returns a named character vector of colors for common tRNA
#' modifications, suitable for use with [plot_tRNA_structure()].
#'
#' @return A named character vector mapping modification names to
#'   hex colors.
#'
#' @export
#'
#' @examples
#' default_mod_palette()
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
