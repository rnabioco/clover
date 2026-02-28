# Generate Sprinzl coordinate files from cmalign output
#
# Aligns tRNA FASTA sequences against the tRNAscan-SE bacterial
# covariance model using cmalign, then assigns Sprinzl positions
# based on the CM consensus structure.
#
# Run with: pixi run Rscript data-raw/sprinzl.R
#
# Prerequisites:
#   - cmalign on PATH
#   - Bacterial CM at data-raw/structures/TRNAinf-bact.cm

devtools::load_all()

script_dir <- "data-raw/structures"
cm_file <- file.path(script_dir, "TRNAinf-bact.cm")
out_dir <- "inst/extdata/sprinzl"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# --- CM-to-Sprinzl mapping table ---
#
# The bacterial tRNAscan-SE CM has 93 consensus columns. Each maps to
# a fixed Sprinzl label and structural region. Insert columns between
# consensus positions receive insertion labels (17a, 20a, 20b, 47).
#
# SS_cons structure (93 consensus + inserts):
# (((((((,,<<<<________>>>>,<<<<<_______>>>>>,,<<<<<<<____>>>>>>>,.,
# <<<<<_______>>>>>)))))))::::

build_cm_mapping <- function() {
  # Returns a data frame with columns: cm_col, sprinzl_label, region
  cm <- data.frame(
    cm_col = 1:93,
    sprinzl_label = NA_character_,
    region = NA_character_,
    stringsAsFactors = FALSE
  )

  assign_range <- function(cm_start, cm_end, sp_start, region) {
    idx <- cm_start:cm_end
    cm$sprinzl_label[idx] <<- as.character(
      sp_start:(sp_start + length(idx) - 1)
    )
    cm$region[idx] <<- region
  }

  assign_range(1, 7, 1, "acceptor-stem")
  cm$region[8:9] <- "unknown"
  cm$sprinzl_label[8:9] <- c("8", "9")
  assign_range(10, 13, 10, "D-stem")
  assign_range(14, 21, 14, "D-loop")
  assign_range(22, 25, 22, "D-stem")
  cm$sprinzl_label[26] <- "26"
  cm$region[26] <- "unknown"
  assign_range(27, 31, 27, "anticodon-stem")
  assign_range(32, 38, 32, "anticodon-loop")
  assign_range(39, 43, 39, "anticodon-stem")
  cm$sprinzl_label[44] <- "44"
  cm$region[44] <- "variable-region"
  cm$sprinzl_label[45] <- "45"
  cm$region[45] <- "variable-region"
  # CM 46-52: variable stem 5' (e-positions, handled per-tRNA)
  cm$region[46:52] <- "variable-arm"
  # CM 53-56: variable loop
  cm$sprinzl_label[53:56] <- paste0("e", 1:4)
  cm$region[53:56] <- "variable-arm"
  # CM 57-63: variable stem 3' (e-positions, handled per-tRNA)
  cm$region[57:63] <- "variable-arm"
  cm$sprinzl_label[64] <- "46"
  cm$region[64] <- "variable-region"
  cm$sprinzl_label[65] <- "48"
  cm$region[65] <- "variable-arm"
  assign_range(66, 70, 49, "T-stem")
  assign_range(71, 77, 54, "T-loop")
  assign_range(78, 82, 61, "T-stem")
  assign_range(83, 89, 66, "acceptor-stem")
  cm$sprinzl_label[90] <- "73"
  cm$region[90] <- "acceptor-tail"
  assign_range(91, 93, 74, "acceptor-tail")

  cm
}

# --- Parse Stockholm alignment ---

parse_stockholm <- function(sto_file) {
  lines <- readLines(sto_file)

  # Extract RF and SS_cons (concatenate across wrapped blocks)
  rf_lines <- grep("^#=GC RF", lines, value = TRUE)
  rf_str <- paste0(sub("^#=GC RF\\s+", "", rf_lines), collapse = "")
  rf_chars <- strsplit(rf_str, "")[[1]]

  ss_lines <- grep("^#=GC SS_cons", lines, value = TRUE)
  ss_str <- paste0(sub("^#=GC SS_cons\\s+", "", ss_lines), collapse = "")
  ss_chars <- strsplit(ss_str, "")[[1]]

  # Extract sequence lines (not comments, not PP, not GC, not blank)
  seq_lines <- lines[!grepl("^#|^//|^$", lines)]
  seq_lines <- seq_lines[!grepl("\\s+PP\\s+", seq_lines)]

  # Concatenate wrapped sequence blocks by name
  seqs <- list()
  for (sl in seq_lines) {
    parts <- strsplit(trimws(sl), "\\s+")[[1]]
    name <- parts[1]
    chunk <- parts[2]
    if (is.null(seqs[[name]])) {
      seqs[[name]] <- chunk
    } else {
      seqs[[name]] <- paste0(seqs[[name]], chunk)
    }
  }

  list(rf = rf_chars, ss = ss_chars, seqs = seqs)
}

# --- Assign Sprinzl positions for one tRNA ---

assign_sprinzl <- function(aln_seq, rf_chars, cm_map) {
  aln_chars <- strsplit(aln_seq, "")[[1]]
  stopifnot(length(aln_chars) == length(rf_chars))

  # First pass: identify CM columns and which have bases
  cm_col <- 0L
  # Track which CM columns in the variable stem have bases
  var_stem5_present <- logical(7) # CM 46-52
  var_stem3_present <- logical(7) # CM 57-63

  for (i in seq_along(aln_chars)) {
    is_cons <- rf_chars[i] != "."
    if (is_cons) {
      cm_col <- cm_col + 1L
    }
    has_base <- toupper(aln_chars[i]) %in% c("A", "C", "G", "U", "T")
    if (has_base && is_cons) {
      if (cm_col >= 46 && cm_col <= 52) {
        var_stem5_present[cm_col - 45] <- TRUE
      }
      if (cm_col >= 57 && cm_col <= 63) {
        var_stem3_present[cm_col - 56] <- TRUE
      }
    }
  }

  # Compute e-position labels for variable stem
  n_stem5 <- sum(var_stem5_present)
  n_stem3 <- sum(var_stem3_present)

  # 5' strand: labels go e11, e12, ..., counting from outer to inner
  # The innermost position (closest to loop) gets the highest index
  # For n bases, labels are: e1(5-n) to e14 if n <= 4
  # For n > 4, extend: e11 to e1n
  stem5_labels <- character(7)
  if (n_stem5 > 0) {
    which_present <- which(var_stem5_present)
    if (n_stem5 <= 4) {
      start_idx <- 5 - n_stem5
    } else {
      start_idx <- 1
    }
    label_idx <- 1
    for (j in which_present) {
      stem5_labels[j] <- paste0("e1", start_idx + label_idx - 1)
      label_idx <- label_idx + 1
    }
  }

  # 3' strand: labels go e2N, e2(N-1), ..., e21 (inner to outer)
  # The innermost position (closest to loop) gets the highest index
  stem3_labels <- character(7)
  if (n_stem3 > 0) {
    which_present <- which(var_stem3_present)
    if (n_stem3 <= 4) {
      start_idx <- 5 - n_stem3
    } else {
      start_idx <- 1
    }
    # 3' strand labels decrease: innermost (leftmost in alignment) gets
    # highest number, outermost (rightmost) gets lowest
    label_idx <- n_stem3
    for (j in which_present) {
      stem3_labels[j] <- paste0("e2", start_idx + label_idx - 1)
      label_idx <- label_idx - 1
    }
  }

  # Second pass: assign Sprinzl positions
  cm_col <- 0L
  seq_pos <- 0L
  prev_cm <- 0L # track previous consensus column for insert labeling
  rows <- list()

  for (i in seq_along(aln_chars)) {
    is_cons <- rf_chars[i] != "."
    if (is_cons) {
      cm_col <- cm_col + 1L
    }
    has_base <- toupper(aln_chars[i]) %in% c("A", "C", "G", "U", "T")

    if (!has_base) {
      if (is_cons) {
        prev_cm <- cm_col
      }
      next
    }

    seq_pos <- seq_pos + 1L
    residue <- toupper(aln_chars[i])
    # Convert U to T for DNA convention
    if (residue == "U") {
      residue <- "T"
    }

    if (is_cons) {
      # Consensus position
      sp_label <- cm_map$sprinzl_label[cm_col]
      sp_region <- cm_map$region[cm_col]

      # Handle variable stem e-positions
      if (cm_col >= 46 && cm_col <= 52 && is.na(sp_label)) {
        sp_label <- stem5_labels[cm_col - 45]
      }
      if (cm_col >= 57 && cm_col <= 63 && is.na(sp_label)) {
        sp_label <- stem3_labels[cm_col - 56]
      }

      rows[[length(rows) + 1]] <- data.frame(
        seq_index = seq_pos,
        sprinzl_label = sp_label,
        global_index = i,
        region = sp_region,
        residue = residue,
        stringsAsFactors = FALSE
      )
      prev_cm <- cm_col
    } else {
      # Insert position: label based on context
      sp_label <- NA_character_
      sp_region <- "unknown"

      if (prev_cm == 17) {
        # Insert between Sprinzl 17 and 18 → 17a
        sp_label <- "17a"
        sp_region <- "D-loop"
      } else if (prev_cm == 20) {
        # Insert between CM 20 and CM 21 → 20a, 20b
        # Check if previous insert in this gap was already 20a
        n_prev_inserts <- sum(
          vapply(
            rows,
            \(r) {
              !is.na(r$sprinzl_label) && r$sprinzl_label %in% c("20a", "20b")
            },
            logical(1)
          )
        )
        if (n_prev_inserts == 0) {
          sp_label <- "20a"
        } else {
          sp_label <- "20b"
        }
        sp_region <- "D-loop"
      } else if (prev_cm == 64) {
        # Insert between CM 64 and CM 65 → position 47
        sp_label <- "47"
        sp_region <- "variable-arm"
      } else {
        # Other inserts: label based on surrounding context
        sp_region <- cm_map$region[min(prev_cm, 93)]
      }

      rows[[length(rows) + 1]] <- data.frame(
        seq_index = seq_pos,
        sprinzl_label = sp_label,
        global_index = i,
        region = sp_region,
        residue = residue,
        stringsAsFactors = FALSE
      )
    }
  }

  dplyr::bind_rows(rows)
}

# --- Run cmalign and generate Sprinzl coords ---

run_cmalign <- function(fasta_path, cm_path) {
  sto_file <- tempfile(fileext = ".sto")
  system2(
    "cmalign",
    args = c("--notrunc", "-o", sto_file, cm_path, fasta_path),
    stdout = "",
    stderr = ""
  )
  sto_file
}

generate_sprinzl_coords <- function(fasta_path, cm_path, trna_prefix = "") {
  # Convert FASTA to RNA if needed
  fa <- Biostrings::readDNAStringSet(fasta_path)
  cm_map <- build_cm_mapping()

  cli::cli_inform(
    "Aligning {length(fa)} sequence{?s} against bacterial CM."
  )

  # Align each sequence individually to avoid batch alignment artifacts
  # (batch alignments can introduce ~ characters in the RF line for
  # elided consensus columns, which breaks CM column counting)
  results <- list()
  for (i in seq_along(fa)) {
    name <- names(fa)[i]
    rna_file <- tempfile(fileext = ".fa")
    rna <- Biostrings::RNAStringSet(fa[i])
    Biostrings::writeXStringSet(rna, rna_file)

    sto_file <- run_cmalign(rna_file, cm_path)
    parsed <- parse_stockholm(sto_file)
    unlink(c(rna_file, sto_file))

    coords <- assign_sprinzl(parsed$seqs[[name]], parsed$rf, cm_map)
    # Build trna_id: strip prefix, convert anticodon DNA→RNA
    trna_id <- sub("^host-|^phage-", "", name)
    trna_id <- sub(
      "(tRNA-[A-Za-z0-9]+-)([ACGT]{3})",
      paste0(
        "\\1",
        chartr(
          "T",
          "U",
          sub(".*tRNA-[A-Za-z0-9]+-([ACGT]{3}).*", "\\1", name)
        )
      ),
      trna_id
    )
    coords$trna_id <- trna_id
    results[[name]] <- coords
  }

  dplyr::bind_rows(results) |>
    dplyr::select(
      trna_id,
      seq_index,
      sprinzl_label,
      global_index,
      region,
      residue
    )
}

# --- Extract T4 phage tRNAs from test data ---

extract_phage_fasta <- function(fasta_path, prefix, output_path) {
  fa <- Biostrings::readDNAStringSet(fasta_path)
  phage <- fa[grep(paste0("^", prefix), names(fa))]
  if (length(phage) == 0) {
    cli::cli_abort("No sequences matching prefix {.val {prefix}} found.")
  }
  Biostrings::writeXStringSet(phage, output_path)
  cli::cli_inform("Wrote {length(phage)} sequence{?s} to {output_path}.")
  output_path
}

# --- Extract T5 phage tRNAs from MODOMICS ---

extract_modomics_fasta <- function(organism, output_path) {
  mod_dict <- load_cached_modifications()
  seqs <- load_cached_sequences(organism)

  if (is.null(seqs) || nrow(seqs) == 0) {
    cli::cli_abort("No MODOMICS sequences for {.val {organism}}.")
  }

  plain_seqs <- vapply(
    seq_len(nrow(seqs)),
    function(i) strip_modifications(seqs$seq[i], mod_dict),
    character(1)
  )

  # Convert RNA (U) to DNA (T) for pipeline consistency
  dna_seqs <- chartr("U", "T", plain_seqs)

  # Build names: phage-tRNA-{AA}-{anticodon_DNA}
  anticodons_dna <- chartr("U", "T", seqs$anticodon)
  # Map 3-letter AA codes to pipeline-style names
  aa_map <- c(
    Ala = "Ala",
    Arg = "Arg",
    Asn = "Asn",
    Asp = "Asp",
    Cys = "Cys",
    Gln = "Gln",
    Glu = "Glu",
    Gly = "Gly",
    His = "His",
    Ile = "Ile",
    Leu = "Leu",
    Lys = "Lys",
    Met = "Met",
    Phe = "Phe",
    Pro = "Pro",
    Ser = "Ser",
    Thr = "Thr",
    Trp = "Trp",
    Tyr = "Tyr",
    Val = "Val",
    Ini = "fMet",
    SeC = "SeC"
  )
  aa_names <- aa_map[seqs$subtype]

  trna_names <- paste0("phage-tRNA-", aa_names, "-", anticodons_dna)

  fa <- Biostrings::DNAStringSet(dna_seqs)
  names(fa) <- trna_names
  Biostrings::writeXStringSet(fa, output_path)
  cli::cli_inform("Wrote {length(fa)} sequence{?s} to {output_path}.")
  output_path
}

# === Generate Sprinzl coordinates ===

fasta_dir <- "data-raw/structures/fasta"
dir.create(fasta_dir, showWarnings = FALSE, recursive = TRUE)

# --- T4 phage tRNAs ---

cli::cli_h1("Phage T4 tRNAs")

t4_fasta <- file.path(fasta_dir, "phageT4-tRNAs.fa")
extract_phage_fasta(
  "inst/extdata/ecoli/trna_only.fa.gz",
  "phage-",
  t4_fasta
)

t4_coords <- generate_sprinzl_coords(t4_fasta, cm_file)
fname <- "phageT4_global_coords.tsv.gz"
readr::write_tsv(t4_coords, file.path(out_dir, fname))
cli::cli_inform(
  "Saved {nrow(t4_coords)} position{?s} for {length(unique(t4_coords$trna_id))} tRNA{?s} to {fname}."
)

# --- T5 phage tRNAs ---

cli::cli_h1("Phage T5 tRNAs")

t5_fasta <- file.path(fasta_dir, "phageT5-tRNAs.fa")
extract_modomics_fasta("Enterobacteria phage T5", t5_fasta)

t5_coords <- generate_sprinzl_coords(t5_fasta, cm_file)
fname <- "phageT5_global_coords.tsv.gz"
readr::write_tsv(t5_coords, file.path(out_dir, fname))
cli::cli_inform(
  "Saved {nrow(t5_coords)} position{?s} for {length(unique(t5_coords$trna_id))} tRNA{?s} to {fname}."
)

# --- E. coli K12 host tRNAs ---

cli::cli_h1("E. coli K12 host tRNAs")

ecoli_fasta <- file.path(fasta_dir, "ecoliK12-tRNAs.fa")
extract_phage_fasta(
  "inst/extdata/ecoli/trna_only.fa.gz",
  "host-",
  ecoli_fasta
)

ecoli_coords <- generate_sprinzl_coords(ecoli_fasta, cm_file)
fname <- "ecoliK12_global_coords.tsv.gz"
readr::write_tsv(ecoli_coords, file.path(out_dir, fname))
cli::cli_inform(
  "Saved {nrow(ecoli_coords)} position{?s} for {length(unique(ecoli_coords$trna_id))} tRNA{?s} to {fname}."
)

cli::cli_h1("Done")
