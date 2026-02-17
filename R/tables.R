# Table functions --------------------------------------------------------------

#' Tabulate top DESeq2 differential expression results.
#'
#' Create a formatted gt table of the top significant tRNAs from
#' [tidy_deseq_results()] output, sorted by p-value.
#'
#' @param data A tibble from [tidy_deseq_results()] with at least
#'   `log2FoldChange`, `pvalue`, `padj`, and `significant` columns.
#' @param lab_col Column name (string) used for row labels. Default
#'   `"ref"`.
#' @param n Maximum number of rows to display. Default `10`.
#'
#' @return An interactive `gt_tbl` object with search, sorting, and
#'   pagination.
#'
#' @export
#'
#' @examples
#' res <- tibble::tibble(
#'   ref = paste0("tRNA-", 1:20),
#'   log2FoldChange = rnorm(20),
#'   pvalue = runif(20, 0, 0.1),
#'   padj = runif(20, 0, 0.2),
#'   significant = c(rep(TRUE, 10), rep(FALSE, 10))
#' )
#' if (requireNamespace("gt", quietly = TRUE)) {
#'   tabulate_deseq(res)
#' }
tabulate_deseq <- function(data, lab_col = "ref", n = 10) {
  rlang::check_installed("gt", reason = "to create formatted tables.")

  tbl_data <- data |>
    dplyr::filter(!is.na(.data$padj)) |>
    dplyr::arrange(.data$pvalue) |>
    utils::head(n) |>
    dplyr::select(
      dplyr::all_of(c(
        lab_col,
        "log2FoldChange",
        "pvalue",
        "padj",
        "significant"
      ))
    )

  tbl_data |>
    gt::gt() |>
    gt::fmt_number(columns = "log2FoldChange", decimals = 2) |>
    gt::fmt_scientific(columns = c("pvalue", "padj")) |>
    gt::cols_label(
      log2FoldChange = "log2 FC",
      pvalue = "p-value",
      padj = "Adjusted p-value"
    ) |>
    gt::opt_row_striping() |>
    gt::opt_interactive(
      use_search = TRUE,
      use_sorting = TRUE,
      use_pagination = TRUE,
      use_page_size_select = TRUE,
      page_size_default = 10
    )
}
