# Identity element functions ---------------------------------------------------

#' Retrieve tRNA aminoacylation identity elements
#'
#' Returns experimentally validated identity elements (determinants
#' and antideterminants) for tRNA aminoacylation, based on data from
#' Giege & Eriani (2023). Elements are defined per amino acid family
#' using Sprinzl numbering.
#'
#' @param organism Character string specifying the organism name
#'   (e.g., `"Escherichia coli"`, `"Saccharomyces cerevisiae"`,
#'   `"Homo sapiens"`). Use [identity_organisms()] to list supported
#'   organisms.
#' @param amino_acid Optional character vector of 3-letter amino acid
#'   codes to filter by (e.g., `"Ala"`, `"Arg"`). If `NULL` (default),
#'   returns elements for all amino acids.
#' @param type Which element types to return: `"both"` (default),
#'   `"determinant"`, or `"antideterminant"`.
#'
#' @return A tibble with columns:
#'   - `amino_acid`: 3-letter amino acid code
#'   - `domain`: evolutionary domain (`"Bacteria"` or `"Eukarya"`)
#'   - `sprinzl_pos`: Sprinzl position number (integer; `NA` for
#'     structural determinants like the long variable arm; `-1` for
#'     the His G\eqn{_{-1}} position)
#'   - `nucleotide`: expected base (A/C/G/U, or `NA` if variable)
#'   - `region`: tRNA structural region name
#'   - `type`: `"determinant"` or `"antideterminant"`
#'   - `strength`: `"strong"` or `"weak"` (determinants only)
#'   - `pair_pos`: partner Sprinzl position for base pairs (`NA`
#'     for single-nucleotide elements)
#'   - `pair_type`: `"wc"` or `"wobble"` (`NA` for unpaired)
#'   - `against_aars`: aaRS blocked (antideterminants only)
#'   - `universal`: `TRUE` if conserved across all domains
#'   - `description`: human-readable description
#'
#' @section Limitations:
#' Selenocysteine tRNA (SeC, anticodon UCA) is not currently
#' supported. Its non-canonical 90-nt structure with an extended
#' variable arm is not represented in the bundled determinant or
#' antideterminant tables, and Sprinzl coordinates for SeC tRNAs are
#' likewise unavailable.
#'
#' @export
#'
#' @references
#' Giege R, Eriani G (2023). "The tRNA identity landscape for
#' aminoacylation and beyond." *Nucleic Acids Research*, 51(4),
#' 1528--1570. \doi{10.1093/nar/gkad007}
#'
#' @examples
#' # All identity elements for E. coli
#' identity_elements("Escherichia coli")
#'
#' # Ala determinants only
#' identity_elements("Escherichia coli", amino_acid = "Ala",
#'   type = "determinant")
identity_elements <- function(
  organism,
  amino_acid = NULL,
  type = c("both", "determinant", "antideterminant")
) {
  type <- match.arg(type)

  org_map <- load_identity_organism_map()
  domain <- org_map$domain[org_map$organism == organism]

  if (length(domain) == 0) {
    supported <- paste0(
      "'",
      org_map$organism,
      "'",
      collapse = ", "
    )
    cli::cli_abort(
      c(
        "Organism {.val {organism}} is not supported.",
        "i" = "Supported organisms: {supported}.",
        "i" = "Use {.fn identity_organisms} to list them."
      )
    )
  }

  result <- dplyr::tibble(
    amino_acid = character(),
    domain = character(),
    sprinzl_pos = integer(),
    nucleotide = character(),
    region = character(),
    type = character(),
    strength = character(),
    pair_pos = integer(),
    pair_type = character(),
    against_aars = character(),
    universal = logical(),
    description = character()
  )

  if (type %in% c("both", "determinant")) {
    dets <- load_cached_determinants()
    dets <- dets[dets$domain == domain, ]
    dets$type <- "determinant"
    dets$against_aars <- NA_character_
    result <- dplyr::bind_rows(result, dets)
  }

  if (type %in% c("both", "antideterminant")) {
    antis <- load_cached_antideterminants()
    antis <- antis[antis$domain == domain, ]
    antis$type <- "antideterminant"
    antis$strength <- NA_character_
    antis$universal <- NA
    result <- dplyr::bind_rows(result, antis)
  }

  result <- result |>
    dplyr::select(
      "amino_acid",
      "domain",
      "sprinzl_pos",
      "nucleotide",
      "region",
      "type",
      "strength",
      "pair_pos",
      "pair_type",
      "against_aars",
      "universal",
      "description"
    )

  if (!is.null(amino_acid)) {
    result <- result[result$amino_acid %in% amino_acid, ]
  }

  result
}

#' List supported organisms for identity elements
#'
#' Returns the names of organisms for which tRNA aminoacylation
#' identity element data is available.
#'
#' @return A character vector of organism names.
#'
#' @export
#'
#' @examples
#' identity_organisms()
identity_organisms <- function() {
  org_map <- load_identity_organism_map()
  org_map$organism
}

#' Retrieve canonical tRNA tertiary contacts
#'
#' Returns the canonical cloverleaf tertiary interactions that
#' stabilize the L-shape and elbow of the tRNA fold: four base
#' triples and five tertiary base pairs. Numbering follows the
#' Sprinzl convention. Modified bases (m7G46, Psi55, m1A58) are
#' encoded as their unmodified equivalents (G, U, A) so callers can
#' compare against in-vitro transcripts directly.
#'
#' Sources: Westhof & Auffinger (2012) and Giege & Eriani (2023).
#' These tertiary contacts are not aaRS identity elements per se but
#' constrain the tRNA fold and may indirectly affect aminoacylation
#' specificity. Use [identity_elements()] for aaRS recognition
#' determinants.
#'
#' @param organism Character string specifying the organism name
#'   (e.g., `"Escherichia coli"`, `"Saccharomyces cerevisiae"`).
#'   Use [identity_organisms()] to list supported organisms.
#'
#' @return A tibble with one row per contact and columns:
#'   - `contact_id`: short label (e.g., `"U8-A14-A21"`, `"G15-C48"`)
#'   - `contact_type`: `"triple"` or `"pair"`
#'   - `domain`: `"Bacteria"` or `"Eukarya"`
#'   - `pos1`, `pos2`, `pos3`: Sprinzl positions (`pos3 = NA` for pairs)
#'   - `nuc1`, `nuc2`, `nuc3`: canonical bases (`nuc3 = NA` for pairs)
#'   - `interaction`: geometry label (e.g., `"reverse-Hoogsteen"`,
#'     `"trans-WC"`, `"Levitt"`)
#'   - `universal`: `TRUE` for contacts conserved across all domains
#'   - `description`: human-readable description
#'
#' @section Limitations:
#' Selenocysteine tRNA (SeC, anticodon UCA) is not currently
#' supported. Its non-canonical 90-nt fold with an extended variable
#' arm does not share the canonical tertiary contact set bundled
#' here, and SeC Sprinzl coordinates are not available.
#'
#' @export
#'
#' @references
#' Westhof E, Auffinger P (2012). "tRNA structure." *Encyclopedia of
#' Life Sciences*. John Wiley & Sons.
#'
#' Giege R, Eriani G (2023). "The tRNA identity landscape for
#' aminoacylation and beyond." *Nucleic Acids Research*, 51(4),
#' 1528--1570. \doi{10.1093/nar/gkad007}
#'
#' @examples
#' tertiary_contacts("Escherichia coli")
#'
#' # Just the four classical triples
#' tc <- tertiary_contacts("Escherichia coli")
#' tc[tc$contact_type == "triple", ]
tertiary_contacts <- function(organism) {
  org_map <- load_identity_organism_map()
  domain <- org_map$domain[org_map$organism == organism]

  if (length(domain) == 0) {
    supported <- paste0(
      "'",
      org_map$organism,
      "'",
      collapse = ", "
    )
    cli::cli_abort(
      c(
        "Organism {.val {organism}} is not supported.",
        "i" = "Supported organisms: {supported}.",
        "i" = "Use {.fn identity_organisms} to list them."
      )
    )
  }

  contacts <- load_cached_tertiary()
  contacts[contacts$domain == domain, ]
}

#' Map identity elements to tRNA sequence positions
#'
#' Converts Sprinzl positions in identity element data to 1-based
#' sequence positions for a specific tRNA, enabling overlay on
#' structure plots via the `outlines` parameter of
#' [plot_tRNA_structure()].
#'
#' @param elements A tibble of identity elements as returned by
#'   [identity_elements()].
#' @param sprinzl_coords A tibble of Sprinzl coordinates as returned
#'   by [read_sprinzl_coords()].
#' @param trna_id Character string identifying the tRNA to map
#'   (must match a `trna_id` value in `sprinzl_coords`).
#'
#' @return The input `elements` tibble with an added `pos` column
#'   containing the 1-based sequence position. Rows where the
#'   Sprinzl position does not exist in the target tRNA (including
#'   His G\eqn{_{-1}} at `sprinzl_pos = -1` and structural
#'   determinants with `sprinzl_pos = NA`) are dropped.
#'
#' @export
#'
#' @examples
#' coords <- read_sprinzl_coords(
#'   clover_example("sprinzl/sacCer_global_coords.tsv.gz")
#' )
#' elems <- identity_elements(
#'   "Saccharomyces cerevisiae",
#'   amino_acid = "Ala"
#' )
#' map_identity_to_trna(elems, coords, "nuc-tRNA-Ala-AGC-1-1")
map_identity_to_trna <- function(elements, sprinzl_coords, trna_id) {
  coords <- sprinzl_coords[sprinzl_coords$trna_id == trna_id, ]

  if (nrow(coords) == 0) {
    cli::cli_abort(
      "tRNA {.val {trna_id}} not found in {.arg sprinzl_coords}."
    )
  }

  # Convert sprinzl_pos to character for joining with sprinzl_label
  elements$sprinzl_int <- as.character(elements$sprinzl_pos)

  merged <- dplyr::inner_join(
    elements,
    coords[, c("sprinzl_label", "pos")],
    by = c("sprinzl_int" = "sprinzl_label")
  )

  merged$sprinzl_int <- NULL
  merged
}

#' Plot tRNA structure with identity element overlays
#'
#' Generates a tRNA cloverleaf structure SVG with aminoacylation
#' identity elements highlighted as colored outlines. Strong
#' determinants are shown in red, weak determinants in blue.
#'
#' @param trna Character string identifying the tRNA
#'   (e.g., `"tRNA-Ala-AGC"`). Use [structure_trnas()] to list
#'   available tRNAs.
#' @param organism Character string specifying the organism name
#'   (e.g., `"Saccharomyces cerevisiae"`).
#' @param sprinzl_coords A tibble of Sprinzl coordinates as returned
#'   by [read_sprinzl_coords()].
#' @param trna_id Character string identifying the tRNA in
#'   `sprinzl_coords` (e.g., `"nuc-tRNA-Ala-AGC-1-1"`). If `NULL`
#'   (default), auto-detected from `trna` by converting the anticodon
#'   to RNA (T to U) and matching the first Sprinzl entry.
#' @param amino_acid Optional 3-letter amino acid code. If `NULL`
#'   (default), extracted from the `trna` name.
#' @param outline_palette Named character vector of colors keyed by
#'   strength (`"strong"`, `"weak"`). Default uses red for strong
#'   and blue for weak determinants.
#' @param ... Additional arguments passed to
#'   [plot_tRNA_structure()].
#'
#' @return The path to the annotated SVG file (invisibly).
#'
#' @export
#'
#' @examples
#' \donttest{
#' coords <- read_sprinzl_coords(
#'   clover_example("sprinzl/sacCer_global_coords.tsv.gz")
#' )
#' plot_identity_structure(
#'   "tRNA-Ala-AGC", "Saccharomyces cerevisiae", coords
#' )
#' }
plot_identity_structure <- function(
  trna,
  organism,
  sprinzl_coords,
  trna_id = NULL,
  amino_acid = NULL,
  outline_palette = NULL,
  ...
) {
  if (is.null(amino_acid)) {
    amino_acid <- strsplit(trna, "-")[[1]][2]
  }

  if (is.null(trna_id)) {
    trna_id <- find_sprinzl_id(trna, sprinzl_coords)
    if (is.null(trna_id)) {
      cli::cli_abort(
        "No Sprinzl coordinate entry found for {.val {trna}}."
      )
    }
  }

  elements <- identity_elements(
    organism,
    amino_acid = amino_acid,
    type = "determinant"
  )

  mapped <- map_identity_to_trna(elements, sprinzl_coords, trna_id)

  if (nrow(mapped) == 0) {
    return(plot_tRNA_structure(trna, organism, ...))
  }

  outlines <- dplyr::tibble(
    pos = mapped$pos,
    group = mapped$strength
  )

  if (is.null(outline_palette)) {
    outline_palette <- c(
      strong = "#E41A1C",
      weak = "#377EB8"
    )
  }

  plot_tRNA_structure(
    trna,
    organism,
    outlines = outlines,
    outline_palette = outline_palette,
    ...
  )
}

#' Plot identity elements for multiple tRNAs in a grid
#'
#' Generates a combined SVG showing tRNA cloverleaf structures
#' arranged in a grid, each annotated with aminoacylation identity
#' elements. Inspired by Figure 2 of Giege & Eriani (2023).
#'
#' @param trnas Character vector of tRNA identifiers
#'   (e.g., `c("tRNA-Ala-AGC", "tRNA-Phe-GAA")`).
#' @param organism Character string specifying the organism name.
#' @param sprinzl_coords A tibble of Sprinzl coordinates as returned
#'   by [read_sprinzl_coords()].
#' @param output Path for the output SVG file. If `NULL` (default),
#'   writes to a temporary file.
#' @param outline_palette Named character vector of colors keyed by
#'   strength. Default uses red for strong and blue for weak.
#' @param ncol Number of columns in the panel grid. If `NULL`
#'   (default), all tRNAs are placed in a single row.
#' @param gap Gap in SVG units between panels. Default 20.
#' @param ... Additional arguments passed to
#'   [plot_identity_structure()].
#'
#' @return The path to the combined SVG file (invisibly).
#'
#' @export
#'
#' @examples
#' \donttest{
#' coords <- read_sprinzl_coords(
#'   clover_example("sprinzl/sacCer_global_coords.tsv.gz")
#' )
#' plot_identity_panel(
#'   c("tRNA-Ala-AGC", "tRNA-Asp-GTC"),
#'   "Saccharomyces cerevisiae",
#'   coords
#' )
#' }
plot_identity_panel <- function(
  trnas,
  organism,
  sprinzl_coords,
  output = NULL,
  outline_palette = NULL,
  ncol = NULL,
  gap = 20,
  ...
) {
  if (is.null(outline_palette)) {
    outline_palette <- c(
      strong = "#E41A1C",
      weak = "#377EB8"
    )
  }

  # Generate individual SVGs
  svg_paths <- vapply(
    trnas,
    function(trna) {
      plot_identity_structure(
        trna,
        organism,
        sprinzl_coords,
        outline_palette = outline_palette,
        position_markers = FALSE,
        ...
      )
    },
    character(1)
  )

  # Read all SVGs and extract dimensions
  svg_docs <- lapply(svg_paths, xml2::read_xml)
  widths <- vapply(
    svg_docs,
    function(doc) {
      vb <- strsplit(
        xml2::xml_attr(xml2::xml_root(doc), "viewBox"),
        "\\s+"
      )[[1]]
      as.numeric(vb[3])
    },
    numeric(1)
  )
  heights <- vapply(
    svg_docs,
    function(doc) {
      vb <- strsplit(
        xml2::xml_attr(xml2::xml_root(doc), "viewBox"),
        "\\s+"
      )[[1]]
      as.numeric(vb[4])
    },
    numeric(1)
  )

  # Grid layout
  n <- length(trnas)
  if (is.null(ncol)) {
    ncol <- n
  }
  nrow <- ceiling(n / ncol)

  # Label height
  label_h <- 16

  # Cell dimensions based on the largest panel
  cell_w <- max(widths)
  cell_h <- max(heights) + label_h

  total_width <- ncol * cell_w + (ncol - 1) * gap
  total_height <- nrow * cell_h + (nrow - 1) * gap

  # Build combined SVG
  combined <- xml2::read_xml(paste0(
    '<svg xmlns="http://www.w3.org/2000/svg" ',
    'viewBox="0 0 ',
    total_width,
    ' ',
    total_height,
    '" ',
    'width="',
    total_width,
    '" height="',
    total_height,
    '"',
    '></svg>'
  ))
  root <- xml2::xml_root(combined)

  for (i in seq_along(svg_docs)) {
    col_idx <- (i - 1) %% ncol
    row_idx <- (i - 1) %/% ncol

    x_offset <- col_idx * (cell_w + gap)
    y_offset <- row_idx * (cell_h + gap)

    panel_g <- xml2::xml_add_child(
      root,
      "g",
      transform = paste0(
        "translate(",
        x_offset,
        ",",
        y_offset + label_h,
        ")"
      )
    )

    # Copy all children from the individual SVG
    children <- xml2::xml_children(xml2::xml_root(svg_docs[[i]]))
    for (child in children) {
      xml2::xml_add_child(panel_g, child)
    }

    # Add amino acid label above each panel
    aa <- strsplit(trnas[i], "-")[[1]][2]
    label_x <- x_offset + widths[i] / 2
    label_node <- xml2::xml_add_child(
      root,
      "text",
      x = as.character(label_x),
      y = as.character(y_offset + label_h - 3),
      "text-anchor" = "middle",
      "font-family" = "Helvetica, Arial, sans-serif",
      "font-size" = "11",
      "font-weight" = "bold",
      fill = "#333333"
    )
    xml2::xml_set_text(label_node, aa)
  }

  # Add shared legend
  legend_g <- xml2::xml_add_child(root, "g")
  legend_y <- total_height - 3

  items <- list()
  for (nm in names(outline_palette)) {
    items[[length(items) + 1]] <- list(
      label = nm,
      color = outline_palette[[nm]]
    )
  }

  legend_x <- total_width / 2 - (length(items) * 60) / 2
  for (item in items) {
    xml2::xml_add_child(
      legend_g,
      "circle",
      cx = as.character(legend_x),
      cy = as.character(legend_y - 3),
      r = "4",
      fill = "none",
      stroke = item$color,
      "stroke-width" = "1.5"
    )
    txt <- xml2::xml_add_child(
      legend_g,
      "text",
      x = as.character(legend_x + 8),
      y = as.character(legend_y),
      "font-family" = "Helvetica, Arial, sans-serif",
      "font-size" = "9",
      fill = "#555555"
    )
    xml2::xml_set_text(txt, item$label)
    legend_x <- legend_x + 60
  }

  if (is.null(output)) {
    output <- tempfile(fileext = ".svg")
  }
  xml2::write_xml(combined, output)
  invisible(output)
}

# Internal helpers -------------------------------------------------------------

find_sprinzl_id <- function(trna, sprinzl_coords) {
  parts <- strsplit(trna, "-")[[1]]
  if (length(parts) >= 3) {
    parts[3] <- chartr("T", "U", parts[3])
  }
  rna_name <- paste(parts, collapse = "-")
  pattern <- paste0("^nuc-", rna_name, "-")

  ids <- unique(sprinzl_coords$trna_id)
  matches <- grep(pattern, ids, value = TRUE)
  if (length(matches) == 0) {
    return(NULL)
  }
  sort(matches)[1]
}

load_cached_determinants <- function() {
  path <- system.file(
    "extdata",
    "identity",
    "determinants.rds",
    package = "clover"
  )
  readRDS(path)
}

load_cached_tertiary <- function() {
  path <- system.file(
    "extdata",
    "identity",
    "tertiary.rds",
    package = "clover"
  )
  readRDS(path)
}

load_cached_antideterminants <- function() {
  path <- system.file(
    "extdata",
    "identity",
    "antideterminants.rds",
    package = "clover"
  )
  readRDS(path)
}

load_identity_organism_map <- function() {
  path <- system.file(
    "extdata",
    "identity",
    "organism_map.rds",
    package = "clover"
  )
  readRDS(path)
}
