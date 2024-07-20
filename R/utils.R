#' Provide working directory for clover example files.
#'
#' @param path path to file
#'
#' @examples
#' clover_example("bcerror.tsv.gz")
#'
#' @export
clover_example <- function(path) {
  system.file("extdata", path, package = "clover", mustWork = TRUE)
}
