#' Map cached MODOMICS modifications onto reference sequences
#'
#' Uses bundled MODOMICS data to map known tRNA modification
#' positions onto user-provided reference sequences using pairwise
#' alignment. No internet connection is required for organisms
#' included in the package (see [modomics_organisms()]). For other
#' organisms, falls back to [fetch_modomics_mods()].
#'
#' @param fasta Path to a FASTA file or a
#'   [Biostrings::DNAStringSet] object containing reference tRNA
#'   sequences.
#' @param organism Character string specifying the organism name
#'   as used in MODOMICS (e.g., `"Saccharomyces cerevisiae"`,
#'   `"Escherichia coli"`).
#' @param min_identity Minimum alignment identity (0--1) required
#'   to accept a match between a MODOMICS sequence and a reference
#'   sequence. Default `0.7`.
#'
#' @return A tibble with columns:
#'   - `ref`: reference sequence name from the FASTA
#'   - `pos`: 1-based position in the reference sequence
#'   - `mod_full`: full modification name
#'     (e.g., "1-methyladenosine")
#'   - `mod1`: short modification name (e.g., "m1A")
#'
#' @section Limitations:
#' Selenocysteine tRNA (SeC, anticodon UCA) is not currently
#' supported. Its non-canonical 90-nt structure with an extended
#' variable arm differs enough from canonical tRNAs that MODOMICS
#' Sec sequences may fail to align at the default `min_identity`
#' threshold; a warning is emitted in that case so the absence is
#' visible rather than silent.
#'
#' @export
#'
#' @examples
#' \donttest{
#' fa <- clover_example("ecoli/trna_only.fa.gz")
#' modomics_mods(fa, "Escherichia coli")
#' }
modomics_mods <- function(fasta, organism, min_identity = 0.7) {
  rlang::check_installed("pwalign")

  if (is.character(fasta)) {
    fasta <- read_fasta(fasta)
  }

  mod_dict <- load_cached_modifications()
  modomics_seqs <- load_cached_sequences(organism)

  if (is.null(modomics_seqs)) {
    cli::cli_inform(
      c(
        "i" = "No cached data for {.val {organism}}.",
        "i" = "Falling back to {.fn fetch_modomics_mods}."
      )
    )
    return(fetch_modomics_mods(
      fasta,
      organism,
      min_identity = min_identity
    ))
  }

  if (nrow(modomics_seqs) == 0) {
    cli::cli_abort(
      "No tRNA sequences found in MODOMICS for {.val {organism}}."
    )
  }

  cli::cli_inform(
    "Processing {nrow(modomics_seqs)} MODOMICS sequence{?s}."
  )

  modomics_entries <- build_modomics_entries(modomics_seqs, mod_dict)

  cli::cli_inform(
    "Matching MODOMICS sequences to reference FASTA."
  )

  result <- match_modomics_to_refs(
    modomics_entries,
    fasta,
    min_identity
  )

  cli::cli_inform(
    "Found {nrow(result)} modification annotation{?s}."
  )
  result
}

#' List organisms with cached MODOMICS data
#'
#' Returns the names of organisms for which MODOMICS tRNA
#' modification data is bundled with the package. These organisms
#' can be used with [modomics_mods()] without internet access.
#'
#' @return A character vector of organism names.
#'
#' @export
#'
#' @examples
#' modomics_organisms()
modomics_organisms <- function() {
  modomics_dir <- system.file(
    "extdata",
    "modomics",
    package = "clover"
  )
  rds_files <- list.files(modomics_dir, pattern = "\\.rds$")
  rds_files <- rds_files[rds_files != "modifications.rds"]
  gsub("_", " ", tools::file_path_sans_ext(rds_files))
}

build_modomics_entries <- function(modomics_seqs, mod_dict) {
  purrr::map(
    seq_len(nrow(modomics_seqs)),
    function(i) {
      mods <- extract_mod_positions(
        modomics_seqs$seq[i],
        mod_dict
      )
      plain_seq <- strip_modifications(
        modomics_seqs$seq[i],
        mod_dict
      )
      list(
        subtype = modomics_seqs$subtype[i],
        anticodon = modomics_seqs$anticodon[i],
        mods = mods,
        plain_seq = plain_seq
      )
    }
  )
}

load_cached_modifications <- function() {
  path <- system.file(
    "extdata",
    "modomics",
    "modifications.rds",
    package = "clover"
  )
  readRDS(path)
}

load_cached_sequences <- function(organism) {
  fname <- paste0(gsub(" ", "_", organism), ".rds")
  path <- system.file(
    "extdata",
    "modomics",
    fname,
    package = "clover"
  )
  if (path == "") {
    return(NULL)
  }
  readRDS(path)
}

#' Fetch tRNA modification annotations from MODOMICS
#'
#' Downloads tRNA modification data from the
#' [MODOMICS](https://genesilico.pl/modomics/) database and maps
#' modification positions onto user-provided reference sequences
#' using pairwise alignment.
#'
#' @param fasta Path to a FASTA file or a
#'   [Biostrings::DNAStringSet] object containing reference tRNA
#'   sequences.
#' @param organism Character string specifying the organism name
#'   as used in MODOMICS (e.g., `"Saccharomyces cerevisiae"`,
#'   `"Escherichia coli"`).
#' @param cache_dir Optional directory path for caching API
#'   responses as RDS files. If `NULL` (default), no caching is
#'   performed.
#' @param min_identity Minimum alignment identity (0--1) required
#'   to accept a match between a MODOMICS sequence and a reference
#'   sequence. Default `0.7`.
#'
#' @return A tibble with columns:
#'   - `ref`: reference sequence name from the FASTA
#'   - `pos`: 1-based position in the reference sequence
#'   - `mod_full`: full modification name
#'     (e.g., "1-methyladenosine")
#'   - `mod1`: short modification name (e.g., "m1A")
#'
#' @export
#'
#' @examples
#' \donttest{
#' fa <- clover_example("ecoli/validated.fa.gz")
#' mods <- fetch_modomics_mods(fa, "Escherichia coli")
#' mods
#' }
fetch_modomics_mods <- function(
  fasta,
  organism,
  cache_dir = NULL,
  min_identity = 0.7
) {
  rlang::check_installed(c("httr2", "jsonlite", "pwalign"))

  if (is.character(fasta)) {
    fasta <- read_fasta(fasta)
  }

  cli::cli_inform("Fetching MODOMICS modification dictionary.")
  mod_dict <- fetch_modomics_modifications(cache_dir)

  cli::cli_inform(
    "Fetching MODOMICS tRNA sequences for {.val {organism}}."
  )
  modomics_seqs <- fetch_modomics_sequences(organism, cache_dir)

  if (nrow(modomics_seqs) == 0) {
    cli::cli_abort(
      "No tRNA sequences found in MODOMICS for {.val {organism}}."
    )
  }

  cli::cli_inform(
    "Processing {nrow(modomics_seqs)} MODOMICS sequence{?s}."
  )

  modomics_entries <- build_modomics_entries(modomics_seqs, mod_dict)

  cli::cli_inform(
    "Matching MODOMICS sequences to reference FASTA."
  )

  result <- match_modomics_to_refs(
    modomics_entries,
    fasta,
    min_identity
  )

  cli::cli_inform(
    "Found {nrow(result)} modification annotation{?s}."
  )
  result
}

fetch_modomics_modifications <- function(cache_dir = NULL) {
  url <- "https://genesilico.pl/modomics/api/modifications?format=json"
  raw_json <- load_or_fetch(url, cache_dir, "modomics_modifications")

  entries <- jsonlite::fromJSON(raw_json, simplifyVector = FALSE)

  tbl <- dplyr::tibble(
    new_abbrev = purrr::map_chr(
      entries,
      \(e) e$new_abbrev %||% NA_character_
    ),
    short_name = purrr::map_chr(
      entries,
      \(e) e$short_name %||% NA_character_
    ),
    name = purrr::map_chr(
      entries,
      \(e) e$name %||% NA_character_
    ),
    reference_moiety = purrr::map_chr(
      entries,
      function(e) {
        rm_val <- e$reference_moiety
        if (is.null(rm_val) || length(rm_val) == 0) {
          return(NA_character_)
        }
        base <- rm_val[[1]]
        if (base %in% c("A", "U", "G", "C", "T")) base else NA_character_
      }
    )
  )

  tbl[!is.na(tbl$new_abbrev) & nchar(tbl$new_abbrev) > 0, ]
}

fetch_modomics_sequences <- function(organism, cache_dir = NULL) {
  url <- paste0(
    "https://genesilico.pl/modomics/api/sequences",
    "?RNAtype=tRNA",
    "&organism=",
    utils::URLencode(organism, reserved = TRUE),
    "&format=json"
  )
  cache_key <- paste0(
    "modomics_sequences_",
    gsub(" ", "_", organism)
  )
  raw_json <- load_or_fetch(url, cache_dir, cache_key)

  entries <- jsonlite::fromJSON(raw_json, simplifyVector = FALSE)

  if (length(entries) == 0) {
    return(dplyr::tibble(
      subtype = character(),
      anticodon = character(),
      seq = character()
    ))
  }

  dplyr::tibble(
    subtype = purrr::map_chr(
      entries,
      \(e) e$subtype %||% NA_character_
    ),
    anticodon = purrr::map_chr(
      entries,
      \(e) e$anticodon %||% NA_character_
    ),
    seq = purrr::map_chr(
      entries,
      \(e) e$seq %||% NA_character_
    )
  )
}

load_or_fetch <- function(url, cache_dir, cache_key) {
  if (!is.null(cache_dir)) {
    cache_file <- file.path(cache_dir, paste0(cache_key, ".rds"))
    if (file.exists(cache_file)) {
      return(readRDS(cache_file))
    }
  }

  resp <- tryCatch(
    httr2::request(url) |>
      httr2::req_user_agent(
        "clover R package (https://github.com/rnabioco/clover)"
      ) |>
      httr2::req_perform(),
    error = function(e) {
      cli::cli_abort(
        c(
          "Failed to fetch data from MODOMICS API.",
          "x" = "URL: {.url {url}}",
          "i" = "Check your internet connection or try again later."
        ),
        parent = e
      )
    }
  )

  data <- httr2::resp_body_string(resp)

  if (!is.null(cache_dir)) {
    dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)
    saveRDS(data, cache_file)
  }

  data
}

extract_mod_positions <- function(seq, mod_dict) {
  if (is.na(seq) || nchar(seq) == 0) {
    return(dplyr::tibble(
      pos = integer(),
      mod_full = character(),
      mod1 = character()
    ))
  }

  chars <- strsplit(seq, "")[[1]]
  standard_bases <- c("A", "U", "G", "C")

  pos <- integer(0)
  mod_full <- character(0)
  mod1 <- character(0)

  for (i in seq_along(chars)) {
    ch <- chars[i]
    if (ch %in% standard_bases) {
      next
    }

    idx <- match(ch, mod_dict$new_abbrev)
    if (!is.na(idx)) {
      pos <- c(pos, i)
      mod_full <- c(mod_full, mod_dict$name[idx])
      mod1 <- c(mod1, mod_dict$short_name[idx])
    }
  }

  dplyr::tibble(pos = pos, mod_full = mod_full, mod1 = mod1)
}

strip_modifications <- function(seq, mod_dict) {
  if (is.na(seq) || nchar(seq) == 0) {
    return("")
  }

  chars <- strsplit(seq, "")[[1]]
  standard_bases <- c("A", "U", "G", "C")

  result <- purrr::map_chr(
    chars,
    function(ch) {
      if (ch %in% standard_bases) {
        return(ch)
      }
      idx <- match(ch, mod_dict$new_abbrev)
      if (!is.na(idx) && !is.na(mod_dict$reference_moiety[idx])) {
        mod_dict$reference_moiety[idx]
      } else {
        "N"
      }
    }
  )

  paste0(result, collapse = "")
}

match_modomics_to_refs <- function(
  modomics_entries,
  fasta,
  min_identity
) {
  fasta_names <- names(fasta)
  results <- vector("list", length(modomics_entries))

  for (i in seq_along(modomics_entries)) {
    entry <- modomics_entries[[i]]

    if (nrow(entry$mods) == 0) {
      next
    }
    if (nchar(entry$plain_seq) == 0) {
      next
    }

    dna_seq <- gsub("U", "T", entry$plain_seq)
    modomics_dna <- Biostrings::DNAString(dna_seq)

    candidate_idx <- find_aa_candidates(
      entry$subtype,
      fasta_names
    )

    if (length(candidate_idx) == 0) {
      next
    }

    candidate_fa <- fasta[candidate_idx]
    n_cand <- length(candidate_fa)

    patterns <- rep(
      Biostrings::DNAStringSet(modomics_dna),
      n_cand
    )
    alns <- pwalign::pairwiseAlignment(
      patterns,
      candidate_fa,
      type = "local"
    )

    scores <- pwalign::score(alns)
    n_matches <- pwalign::nmatch(alns)
    aln_widths <- nchar(as.character(
      pwalign::alignedPattern(alns)
    ))
    identities <- n_matches / aln_widths

    valid <- identities >= min_identity
    if (!any(valid)) {
      cli::cli_warn(
        c(
          "!" = paste0(
            "No reference matched MODOMICS ",
            "{.val {entry$subtype}-{entry$anticodon}} ",
            "(length {nchar(entry$plain_seq)}) at ",
            "min_identity = {min_identity}; best identity was ",
            "{round(max(identities), 2)}."
          )
        )
      )
      next
    }

    valid_scores <- ifelse(valid, scores, -Inf)
    best <- which.max(valid_scores)

    best_aln <- alns[best]
    best_ref <- fasta_names[candidate_idx[best]]

    transferred <- transfer_positions(best_aln, entry$mods)
    if (nrow(transferred) > 0) {
      transferred$ref <- best_ref
      results[[i]] <- transferred
    }
  }

  result <- dplyr::bind_rows(results)

  if (nrow(result) == 0) {
    return(dplyr::tibble(
      ref = character(),
      pos = integer(),
      mod_full = character(),
      mod1 = character()
    ))
  }

  result |>
    dplyr::select(ref, pos, mod_full, mod1)
}

find_aa_candidates <- function(subtype, fasta_names) {
  patterns <- subtype
  if (subtype == "Ini") {
    patterns <- c("Ini", "iMet", "Met")
  }

  candidate_idx <- integer(0)
  for (pat in patterns) {
    aa_regex <- paste0("[-_]", pat, "\\d*[-_(]")
    idx <- grep(aa_regex, fasta_names, ignore.case = TRUE)
    candidate_idx <- union(candidate_idx, idx)
  }

  candidate_idx
}

transfer_positions <- function(alignment, mod_positions) {
  aln_pattern_str <- as.character(
    pwalign::alignedPattern(alignment)
  )
  aln_subject_str <- as.character(
    pwalign::alignedSubject(alignment)
  )

  p_chars <- strsplit(aln_pattern_str, "")[[1]]
  s_chars <- strsplit(aln_subject_str, "")[[1]]

  p_start <- Biostrings::start(pwalign::pattern(alignment))
  s_start <- Biostrings::start(pwalign::subject(alignment))

  p_pos <- p_start - 1L
  s_pos <- s_start - 1L

  pos_map <- integer(0)
  names_vec <- integer(0)

  for (k in seq_along(p_chars)) {
    p_is_base <- p_chars[k] != "-"
    s_is_base <- s_chars[k] != "-"

    if (p_is_base) {
      p_pos <- p_pos + 1L
    }
    if (s_is_base) {
      s_pos <- s_pos + 1L
    }

    if (p_is_base && s_is_base) {
      pos_map <- c(pos_map, s_pos)
      names_vec <- c(names_vec, p_pos)
    }
  }

  names(pos_map) <- as.character(names_vec)

  mapped_pos <- pos_map[as.character(mod_positions$pos)]

  result <- mod_positions
  result$pos <- as.integer(unname(mapped_pos))

  result[!is.na(result$pos), ]
}
