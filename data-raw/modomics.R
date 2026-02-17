# Fetch and cache MODOMICS data for bundled organisms
#
# Run with: Rscript data-raw/modomics.R

devtools::load_all()

organisms <- c(
  "Escherichia coli",
  "Saccharomyces cerevisiae",
  "Homo sapiens",
  "Mus musculus",
  "Drosophila melanogaster",
  "Caenorhabditis elegans"
)

out_dir <- "inst/extdata/modomics"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Cache modification dictionary
cli::cli_inform("Fetching MODOMICS modification dictionary.")
mod_dict <- fetch_modomics_modifications()
saveRDS(mod_dict, file.path(out_dir, "modifications.rds"))
cli::cli_inform(
  "Saved {nrow(mod_dict)} modification{?s} to modifications.rds."
)

# Cache per-organism tRNA sequences
for (org in organisms) {
  cli::cli_inform("Fetching tRNA sequences for {.val {org}}.")
  seqs <- fetch_modomics_sequences(org)

  if (nrow(seqs) == 0) {
    cli::cli_warn("No tRNA sequences found for {.val {org}}, skipping.")
    next
  }

  fname <- paste0(gsub(" ", "_", org), ".rds")
  saveRDS(seqs, file.path(out_dir, fname))
  cli::cli_inform("Saved {nrow(seqs)} sequence{?s} to {fname}.")
}

cli::cli_inform("Done.")
