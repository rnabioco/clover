# PCoA and rewiring analysis functions ------------------------------------------

#' Prepare a matrix for rewiring PCoA analysis.
#'
#' Build a wide matrix from isodecoder-level relative odds ratios
#' suitable for PCoA or distance calculations. Rows are isodecoders,
#' columns are position pair comparisons.
#'
#' @param data A tibble with columns `isodecoder`, `pos1`, `pos2`, and
#'   `ror`. Typically the output of [compute_ror_isodecoder()].
#' @param sig_only Logical; filter to rows where `significant` is
#'   `TRUE`? Default `TRUE`.
#' @param value_cap Maximum absolute ROR value. Values beyond this are
#'   capped. Default `10`.
#'
#' @return A numeric matrix with isodecoders as row names and position
#'   pair labels (`pos1_vs_pos2`) as column names.
#'
#' @export
#'
#' @examples
#' df <- tibble::tibble(
#'   isodecoder = rep(c("tRNA-Ala", "tRNA-Gly"), each = 2),
#'   pos1 = c(20, 34, 20, 34),
#'   pos2 = c(34, 58, 34, 58),
#'   ror = c(1.5, -0.8, 0.3, 2.1),
#'   significant = c(TRUE, TRUE, TRUE, TRUE)
#' )
#' prepare_rewiring_matrix(df)
prepare_rewiring_matrix <- function(data, sig_only = TRUE, value_cap = 10) {
  if (sig_only) {
    data <- dplyr::filter(data, significant)
  }

  data <- dplyr::mutate(
    data,
    comparison = paste0(pos1, "_vs_", pos2)
  )

  # Aggregate duplicates
  data_agg <- data |>
    dplyr::group_by(isodecoder, comparison) |>
    dplyr::summarise(ror = mean(ror, na.rm = TRUE), .groups = "drop")

  mat <- data_agg |>
    tidyr::pivot_wider(
      names_from = comparison,
      values_from = ror,
      values_fill = 0
    ) |>
    tibble::column_to_rownames("isodecoder") |>
    as.matrix()

  mat[mat > value_cap] <- value_cap
  mat[mat < -value_cap] <- -value_cap

  mat
}

#' Calculate rewiring scores from a ROR matrix.
#'
#' Summarize the magnitude and extent of rewiring for each isodecoder
#' using Euclidean magnitude, mean and max absolute change, and the
#' number of non-zero comparisons.
#'
#' @param mat A numeric matrix from [prepare_rewiring_matrix()], with
#'   isodecoders as row names.
#'
#' @return A tibble with columns: `isodecoder`,
#'   `euclidean_magnitude`, `mean_abs_change`, `max_abs_change`, and
#'   `n_nonzero`, sorted by descending `euclidean_magnitude`.
#'
#' @export
#'
#' @examples
#' mat <- matrix(
#'   c(1.5, -0.8, 0.3, 2.1),
#'   nrow = 2,
#'   dimnames = list(c("tRNA-Ala", "tRNA-Gly"), c("20_vs_34", "34_vs_58"))
#' )
#' calculate_rewiring_scores(mat)
calculate_rewiring_scores <- function(mat) {
  tibble::tibble(
    isodecoder = rownames(mat),
    euclidean_magnitude = sqrt(rowSums(mat^2)),
    mean_abs_change = rowMeans(abs(mat)),
    max_abs_change = apply(abs(mat), 1, max),
    n_nonzero = rowSums(mat != 0)
  ) |>
    dplyr::arrange(dplyr::desc(euclidean_magnitude))
}

#' Perform PCoA on a rewiring matrix.
#'
#' Run classical multidimensional scaling (PCoA) on a Euclidean
#' distance matrix derived from the ROR matrix.
#'
#' @param mat A numeric matrix from [prepare_rewiring_matrix()].
#' @param k Number of dimensions. Default `2`.
#'
#' @return A list with elements:
#'   - `coordinates`: a tibble with `isodecoder`, `PC1`, `PC2`, ...
#'   - `variance_explained`: numeric vector of percent variance
#'     explained for each dimension
#'   - `eigenvalues`: full eigenvalue vector from [stats::cmdscale()]
#'
#' @export
#'
#' @examples
#' mat <- matrix(
#'   c(1.5, -0.8, 0.3, 2.1, 0.5, -1.2),
#'   nrow = 3,
#'   dimnames = list(
#'     c("tRNA-Ala", "tRNA-Gly", "tRNA-Ser"),
#'     c("20_vs_34", "34_vs_58")
#'   )
#' )
#' perform_pcoa(mat)
perform_pcoa <- function(mat, k = 2) {
  dist_mat <- stats::dist(mat, method = "euclidean")
  pcoa <- stats::cmdscale(dist_mat, k = k, eig = TRUE)

  coords <- as.data.frame(pcoa$points)
  colnames(coords) <- paste0("PC", seq_len(k))
  coords$isodecoder <- rownames(mat)

  coords <- tibble::as_tibble(coords)

  # Calculate variance explained
  eig <- pcoa$eig
  var_explained <- eig / sum(abs(eig)) * 100

  list(
    coordinates = coords,
    variance_explained = var_explained[seq_len(k)],
    eigenvalues = eig
  )
}
