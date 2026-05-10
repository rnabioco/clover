# Encode canonical tRNA tertiary contacts (base triples and tertiary
# base pairs) that stabilize the L-shape and elbow of the tRNA fold.
#
# Sources:
#   Westhof E, Auffinger P (2012) "tRNA structure" Encyclopedia of Life
#     Sciences, John Wiley & Sons.
#   Giege R, Eriani G (2023) "The tRNA identity landscape for
#     aminoacylation and beyond" Nucleic Acids Research, gkad007 (Fig 1B).
#   Numbering follows Sprinzl convention.
#
# Modified bases are reduced to their unmodified equivalents for the
# canonical-base comparison: m7G46 -> G, Psi55 -> U, m1A58 -> A.
#
# Run with: Rscript data-raw/tertiary.R

library(dplyr, warn.conflicts = FALSE)

out_dir <- "inst/extdata/identity"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# -- Helpers --
triple <- function(
  contact_id,
  domain,
  pos1,
  pos2,
  pos3,
  nuc1,
  nuc2,
  nuc3,
  interaction,
  universal = FALSE,
  description = NA_character_
) {
  tibble(
    contact_id = contact_id,
    contact_type = "triple",
    domain = domain,
    pos1 = as.integer(pos1),
    pos2 = as.integer(pos2),
    pos3 = as.integer(pos3),
    nuc1 = nuc1,
    nuc2 = nuc2,
    nuc3 = nuc3,
    interaction = interaction,
    universal = universal,
    description = description
  )
}

tpair <- function(
  contact_id,
  domain,
  pos1,
  pos2,
  nuc1,
  nuc2,
  interaction,
  universal = FALSE,
  description = NA_character_
) {
  tibble(
    contact_id = contact_id,
    contact_type = "pair",
    domain = domain,
    pos1 = as.integer(pos1),
    pos2 = as.integer(pos2),
    pos3 = NA_integer_,
    nuc1 = nuc1,
    nuc2 = nuc2,
    nuc3 = NA_character_,
    interaction = interaction,
    universal = universal,
    description = description
  )
}

# Build the same 9 contacts for both domains. Bacterial and eukaryal
# canonical tRNAs share the cloverleaf tertiary core; deviations
# (e.g. type-II tRNAs with long variable arms) are captured downstream
# by per-isoacceptor intactness checks rather than separate domain rows.
build_for_domain <- function(d) {
  bind_rows(
    # --- Triples ---
    triple(
      "U8-A14-A21",
      d,
      8,
      14,
      21,
      "U",
      "A",
      "A",
      "reverse-Hoogsteen",
      universal = TRUE,
      description = "U8-A14-A21 elbow triple (reverse-Hoogsteen + Watson-Crick)"
    ),
    triple(
      "A9-U12-A23",
      d,
      9,
      12,
      23,
      "A",
      "U",
      "A",
      "trans-WC",
      universal = TRUE,
      description = "9-12-23 augmented D-stem triple"
    ),
    triple(
      "G10-C25-G45",
      d,
      10,
      25,
      45,
      "G",
      "C",
      "G",
      "trans-WC",
      universal = TRUE,
      description = "G10-C25-G45 D-stem/variable-loop triple"
    ),
    triple(
      "C13-G22-G46",
      d,
      13,
      22,
      46,
      "C",
      "G",
      "G",
      "trans-WC",
      universal = TRUE,
      description = "C13-G22-m7G46 D-stem/variable-loop triple"
    ),

    # --- Tertiary pairs ---
    tpair(
      "G15-C48",
      d,
      15,
      48,
      "G",
      "C",
      "Levitt",
      description = "G15-C48 Levitt pair (D-loop / variable-loop)"
    ),
    tpair(
      "G18-U55",
      d,
      18,
      55,
      "G",
      "U",
      "trans-WC",
      description = "G18-Psi55 D-loop / T-loop tertiary pair (Psi55 reads as U)"
    ),
    tpair(
      "G19-C56",
      d,
      19,
      56,
      "G",
      "C",
      "WC",
      description = "G19-C56 D-loop / T-loop tertiary pair"
    ),
    tpair(
      "A26-G44",
      d,
      26,
      44,
      "A",
      "G",
      "non-canonical",
      description = "A26-G44 D-stem / variable-loop tertiary pair"
    ),
    tpair(
      "U54-A58",
      d,
      54,
      58,
      "U",
      "A",
      "reverse-Hoogsteen",
      universal = TRUE,
      description = "T54-m1A58 T-loop reverse-Hoogsteen pair (T = U; m1A reads as A)"
    )
  )
}

tertiary <- bind_rows(
  build_for_domain("Bacteria"),
  build_for_domain("Eukarya")
)

cli::cli_inform("Saving tertiary contacts ({nrow(tertiary)} rows).")
saveRDS(tertiary, file.path(out_dir, "tertiary.rds"))

cli::cli_inform("Done.")
