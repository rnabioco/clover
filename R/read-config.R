# Pipeline config reading functions -------------------------------------------

#' Read a Snakemake pipeline configuration file.
#'
#' Parse a YAML config file from the tRNA sequencing pipeline and return
#' a structured list with sample information, output directory, and FASTA
#' reference path. Relative paths are resolved relative to the config file
#' directory.
#'
#' @param config_path Path to a Snakemake `config.yaml` file.
#'
#' @return A list with elements:
#'   - `samples`: A tibble with columns `sample_id` and `data_path`.
#'   - `output_dir`: Resolved path to pipeline output directory.
#'   - `fasta`: Resolved path to reference FASTA file.
#'   - `config_dir`: Directory containing the config file.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' config <- read_pipeline_config("path/to/config.yaml")
#' config$samples
#' config$output_dir
#' }
read_pipeline_config <- function(config_path) {
  config_path <- normalizePath(config_path, mustWork = TRUE)
  config_dir <- dirname(config_path)

  cfg <- yaml::read_yaml(config_path)

  # Parse samples: either inline YAML list or external TSV file
  samples <- .parse_samples(cfg, config_dir)

  # Resolve paths relative to config directory
  # Support both output_dir and output_directory keys
  output_dir_raw <- cfg[["output_dir"]] %||% cfg[["output_directory"]]
  if (is.null(output_dir_raw)) {
    stop("Config file must contain 'output_dir' or 'output_directory'.",
      call. = FALSE
    )
  }
  output_dir <- .resolve_path(output_dir_raw, config_dir)
  fasta <- .resolve_path(cfg[["fasta"]], config_dir)

  list(
    samples = samples,
    output_dir = output_dir,
    fasta = fasta,
    config_dir = config_dir
  )
}

#' List pipeline output files for each sample.
#'
#' Given a parsed pipeline configuration, construct expected file paths for
#' each sample and result type.
#'
#' @param config A list returned by [read_pipeline_config()].
#' @param types Character vector of result types. Valid values:
#'   `"charging"`, `"bcerror"`, `"odds_ratios"`, `"align_stats"`.
#'
#' @return A tibble with columns `sample_id`, `type`, and `path`.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' config <- read_pipeline_config("path/to/config.yaml")
#' list_pipeline_files(config, types = c("charging", "bcerror"))
#' }
list_pipeline_files <- function(
  config,
  types = c("charging", "bcerror", "odds_ratios", "align_stats")
) {
  types <- match.arg(types, several.ok = TRUE)

  suffix_map <- c(
    charging = "charging.cpm",
    bcerror = "bcerror",
    odds_ratios = "odds_ratios",
    align_stats = "align_stats"
  )

  sample_ids <- config$samples$sample_id
  output_dir <- config$output_dir

  rows <- lapply(types, function(type) {
    suffix <- suffix_map[[type]]
    tibble::tibble(
      sample_id = sample_ids,
      type = type,
      path = file.path(
        output_dir,
        "summary",
        "tables",
        sample_ids,
        paste0(sample_ids, ".", suffix, ".tsv.gz")
      )
    )
  })

  dplyr::bind_rows(rows)
}

#' Read pipeline results for all samples.
#'
#' Convenience function that reads a pipeline config, discovers output files,
#' and loads data for the requested result types. Each result type is returned
#' as a combined tibble with a `sample_id` column.
#'
#' @param config_path Path to a Snakemake `config.yaml` file.
#' @param types Character vector of result types to load. Valid values:
#'   `"charging"`, `"bcerror"`, `"odds_ratios"`.
#'
#' @return A named list of tibbles, one per requested type. Each tibble
#'   includes a `sample_id` column identifying the source sample.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' results <- read_pipeline_results("path/to/config.yaml")
#' results$charging
#' results$odds_ratios
#' }
read_pipeline_results <- function(
  config_path,
  types = c("charging", "bcerror", "odds_ratios")
) {
  types <- match.arg(types, several.ok = TRUE)

  config <- read_pipeline_config(config_path)
  files <- list_pipeline_files(config, types = types)

  reader_map <- list(
    charging = read_charging,
    bcerror = read_bcerror,
    odds_ratios = read_odds_ratios
  )

  results <- lapply(types, function(type) {
    type_files <- files[files$type == type, ]
    paths <- stats::setNames(type_files$path, type_files$sample_id)

    # Skip types where no files exist
    existing <- vapply(paths, file.exists, logical(1))
    if (!any(existing)) {
      return(NULL)
    }
    paths <- paths[existing]

    if (type == "charging") {
      read_charging_multi(paths)
    } else if (type == "odds_ratios") {
      read_odds_ratios_multi(paths)
    } else {
      # Generic reader: read each file and bind with sample_id
      tbls <- lapply(names(paths), function(sid) {
        tbl <- reader_map[[type]](paths[[sid]])
        tbl$sample_id <- sid
        tbl
      })
      dplyr::bind_rows(tbls)
    }
  })

  stats::setNames(results, types)
}

# Internal helpers -----------------------------------------------------------

#' Read a samples TSV file, handling both headered and headerless formats.
#'
#' The pipeline's samples.tsv is headerless (two columns: sample_id, data_path).
#' This function auto-detects whether the file has a header by checking if
#' the first line contains "sample_id".
#' @noRd
.read_samples_tsv <- function(path) {
  first_line <- readLines(path, n = 1)
  has_header <- grepl("sample_id", first_line, fixed = TRUE)

  if (has_header) {
    readr::read_tsv(path, show_col_types = FALSE)
  } else {
    readr::read_tsv(
      path,
      col_names = c("sample_id", "data_path"),
      show_col_types = FALSE
    )
  }
}

#' Resolve a path relative to a base directory.
#' @noRd
.resolve_path <- function(path, base_dir) {
  if (is.null(path)) {
    return(NULL)
  }
  path <- as.character(path)
  if (startsWith(path, "/") || startsWith(path, "~")) {
    normalizePath(path, mustWork = FALSE)
  } else {
    normalizePath(file.path(base_dir, path), mustWork = FALSE)
  }
}

#' Parse sample information from config.
#' @noRd
.parse_samples <- function(cfg, config_dir) {
  samples_entry <- cfg[["samples"]]

  if (is.null(samples_entry)) {
    stop("Config file must contain a 'samples' entry.", call. = FALSE)
  }

  # If samples is a string, treat as path to a TSV file
  if (is.character(samples_entry) && length(samples_entry) == 1) {
    samples_path <- .resolve_path(samples_entry, config_dir)
    return(.read_samples_tsv(samples_path))
  }

  # If samples is a list/map, convert to tibble
  if (is.list(samples_entry)) {
    tibble::tibble(
      sample_id = names(samples_entry),
      data_path = vapply(
        samples_entry,
        function(x) .resolve_path(as.character(x), config_dir),
        character(1)
      )
    )
  } else {
    stop("Unsupported 'samples' format in config file.", call. = FALSE)
  }
}
