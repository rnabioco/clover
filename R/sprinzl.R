# Sprinzl coordinate functions ------------------------------------------------

#' Read Sprinzl coordinates from a global coordinates TSV file.
#'
#' Read a tRNAs-in-space global coordinates file and return a tibble with
#' standardized column names. These files map each position of each tRNA
#' to its Sprinzl numbering, structural region, and global alignment index.
#'
#' @param path Path to a global coordinates TSV file (may be gzipped).
#'
#' @return A tibble with columns:
#'   - `trna_id`: tRNA identifier (e.g., `nuc-tRNA-Ala-AGC-1-1`)
#'   - `pos`: 1-based position in tRNA body
#'   - `sprinzl_label`: Sprinzl position (character: "1", "20a", "47:e1", etc.)
#'   - `global_index`: cross-tRNA universal alignment position
#'   - `region`: structural region (acceptor-stem, D-stem, D-loop, etc.)
#'   - `residue`: reference nucleotide
#'
#' @export
#'
#' @examples
#' path <- clover_example("sprinzl/sacCer_global_coords.tsv.gz")
#' coords <- read_sprinzl_coords(path)
#' coords
read_sprinzl_coords <- function(path) {
  readr::read_tsv(
    path,
    col_types = readr::cols(sprinzl_label = readr::col_character()),
    show_col_types = FALSE
  ) |>
    dplyr::rename("pos" = "seq_index") |>
    dplyr::select(
      "trna_id",
      "pos",
      "sprinzl_label",
      "global_index",
      "region",
      "residue"
    )
}

#' Order Sprinzl position labels.
#'
#' Create an ordered factor from Sprinzl position labels. Handles numeric
#' positions and lettered insertions (e.g., "20a", "20A", "47:e1").
#'
#' @param labels Character vector of Sprinzl labels.
#'
#' @return A factor with levels ordered by position number then suffix.
#'
#' @export
#'
#' @examples
#' order_sprinzl_positions(c("20a", "20", "1", "21"))
order_sprinzl_positions <- function(labels) {
  ordering <- tibble::tibble(label = unique(labels)) |>
    dplyr::filter(!is.na(.data$label), .data$label != "") |>
    dplyr::mutate(
      # Variable arm extra positions (e.g., e11, e1, e21) go between 45-46
      .is_var_extra = stringr::str_detect(.data$label, "^e\\d+$"),
      base_num = dplyr::case_when(
        .data$.is_var_extra ~ 45,
        .default = as.numeric(stringr::str_extract(.data$label, "^\\d+"))
      ),
      suffix = stringr::str_extract(.data$label, "[a-zA-Z]$|:e\\d+$"),
      suffix_order = dplyr::case_when(
        .data$.is_var_extra ~
          500 + as.numeric(stringr::str_extract(.data$label, "\\d+$")),
        is.na(.data$suffix) ~ 0,
        stringr::str_detect(.data$suffix, "^[a-z]$") ~
          match(.data$suffix, letters),
        stringr::str_detect(.data$suffix, "^[A-Z]$") ~
          match(.data$suffix, LETTERS),
        stringr::str_detect(.data$suffix, "^:e") ~
          100 + as.numeric(stringr::str_extract(.data$suffix, "\\d+$")),
        .default = 0
      ),
      sort_key = .data$base_num + .data$suffix_order / 1000
    ) |>
    dplyr::arrange(.data$sort_key)

  factor(labels, levels = ordering$label)
}

#' tRNA structural regions.
#'
#' Return a named list mapping canonical tRNA structural region names
#' to integer vectors of Sprinzl position numbers.
#'
#' @return A named list of integer vectors.
#'
#' @export
#'
#' @examples
#' trna_regions()
trna_regions <- function() {
  list(
    acceptor_stem = c(1L:7L, 66L:72L),
    d_arm = 8L:25L,
    anticodon_stem = c(26L:31L, 39L:44L),
    anticodon_loop = 32L:38L,
    variable_loop = 44L:48L,
    t_arm = 49L:65L,
    discriminator = 73L,
    cca = 74L:76L
  )
}

#' Compare a query tRNA to a reference by Sprinzl label.
#'
#' Aligns two tRNAs through their Sprinzl position labels and classifies
#' each reference position as `"match"`, `"substitution"` (same label,
#' different base) or `"deletion"` (label absent from the query). Query
#' positions whose label is absent from the reference are returned
#' separately as insertions.
#'
#' @param reference_coords,query_coords Tibbles of Sprinzl coordinates as
#'   returned by [read_sprinzl_coords()].
#' @param reference_id,query_id The `trna_id` of the reference and query tRNA
#'   within their coordinate tables.
#'
#' @return A list with two tibbles:
#'   - `reference`: one row per reference position with columns `pos`,
#'     `sprinzl_label`, `ref_base`, `query_pos`, `query_base` and `status`.
#'   - `insertions`: one row per query-only position with columns
#'     `query_pos`, `sprinzl_label`, `query_base`, `after_pos` and
#'     `after_label`, the reference position (sequence index and Sprinzl
#'     label) the inserted base follows. Both are `0`/`NA` at the 5' end.
#'
#' @export
#'
#' @examples
#' ecoli <- read_sprinzl_coords(
#'   clover_example("sprinzl/ecoliK12_global_coords.tsv.gz")
#' )
#' t4 <- read_sprinzl_coords(
#'   clover_example("sprinzl/phageT4_global_coords.tsv.gz")
#' )
#' d <- diff_sprinzl(ecoli, "tRNA-Gly-UCC-1-1", t4, "tRNA-Gly-UCC")
#' d$reference[d$reference$status != "match", ]
diff_sprinzl <- function(
  reference_coords,
  reference_id,
  query_coords,
  query_id
) {
  ref <- sprinzl_trna(reference_coords, reference_id, "reference_coords")
  qry <- sprinzl_trna(query_coords, query_id, "query_coords")

  idx <- match(ref$sprinzl_label, qry$sprinzl_label)
  status <- dplyr::case_when(
    is.na(idx) ~ "deletion",
    qry$base[idx] == ref$base ~ "match",
    .default = "substitution"
  )
  reference <- tibble::tibble(
    pos = ref$pos,
    sprinzl_label = ref$sprinzl_label,
    ref_base = ref$base,
    query_pos = qry$pos[idx],
    query_base = qry$base[idx],
    status = status
  )

  # The reference position each query row follows: the last query row at
  # or before it that also exists in the reference (0 at the 5' end).
  ref_pos_of_query <- ref$pos[match(qry$sprinzl_label, ref$sprinzl_label)]
  after <- cummax(replace(ref_pos_of_query, is.na(ref_pos_of_query), 0L))
  only_query <- is.na(ref_pos_of_query)
  insertions <- tibble::tibble(
    query_pos = qry$pos[only_query],
    sprinzl_label = qry$sprinzl_label[only_query],
    query_base = qry$base[only_query],
    after_pos = as.integer(after[only_query]),
    after_label = ref$sprinzl_label[match(after[only_query], ref$pos)]
  )

  list(reference = reference, insertions = insertions)
}

sprinzl_trna <- function(coords, trna_id, arg) {
  rows <- coords[coords$trna_id == trna_id, ]
  if (nrow(rows) == 0) {
    cli::cli_abort("tRNA {.val {trna_id}} not found in {.arg {arg}}.")
  }
  if (anyDuplicated(rows$sprinzl_label)) {
    cli::cli_abort(
      "{.val {trna_id}} has duplicated Sprinzl labels in {.arg {arg}}."
    )
  }
  rows <- rows[order(rows$pos), ]
  tibble::tibble(
    pos = as.integer(rows$pos),
    sprinzl_label = rows$sprinzl_label,
    base = chartr("T", "U", toupper(rows$residue))
  )
}
