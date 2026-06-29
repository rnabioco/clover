# Encode tRNA aminoacylation identity elements from Giege & Eriani (2023)
#
# Source: Giege R, Eriani G (2023) "The tRNA identity landscape for
# aminoacylation and beyond" Nucleic Acids Research, gkad007.
# Tables 1, 3, 4 and Figure 2.
#
# Run with: Rscript data-raw/identity.R

library(dplyr, warn.conflicts = FALSE)

out_dir <- "inst/extdata/identity"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# -- Helper to build determinant rows --
det <- function(
  aa,
  domain,
  sprinzl_pos,
  nucleotide = NA_character_,
  region = NA_character_,
  strength = "strong",
  pair_pos = NA_integer_,
  pair_type = NA_character_,
  universal = FALSE,
  description = NA_character_
) {
  tibble(
    amino_acid = aa,
    domain = domain,
    sprinzl_pos = as.integer(sprinzl_pos),
    nucleotide = nucleotide,
    region = region,
    strength = strength,
    pair_pos = as.integer(pair_pos),
    pair_type = pair_type,
    universal = universal,
    description = description
  )
}

# -- Helper to build antideterminant rows --
anti <- function(
  aa,
  domain,
  sprinzl_pos,
  nucleotide = NA_character_,
  region = NA_character_,
  against_aars = NA_character_,
  pair_pos = NA_integer_,
  pair_type = NA_character_,
  description = NA_character_
) {
  tibble(
    amino_acid = aa,
    domain = domain,
    sprinzl_pos = as.integer(sprinzl_pos),
    nucleotide = nucleotide,
    region = region,
    against_aars = against_aars,
    pair_pos = as.integer(pair_pos),
    pair_type = pair_type,
    description = description
  )
}

# -- Determinants (Tables 1 and 4, Figure 2) --
# Encoding all 20 amino acid families x 2 domains (Bacteria, Eukarya)
# Strength: "strong" = major determinant, "weak" = minor determinant
# Universal (Table 4): TRUE if conserved across all three domains

determinants <- bind_rows(
  # === Ala ===
  det(
    "Ala",
    "Bacteria",
    73,
    "A",
    "discriminator",
    "strong",
    description = "A73 discriminator"
  ),

  det(
    "Ala",
    "Bacteria",
    3,
    "G",
    "acceptor_stem",
    "strong",
    pair_pos = 70L,
    pair_type = "wobble",
    universal = TRUE,
    description = "G3-U70 wobble pair"
  ),
  det(
    "Ala",
    "Bacteria",
    70,
    "U",
    "acceptor_stem",
    "strong",
    pair_pos = 3L,
    pair_type = "wobble",
    universal = TRUE,
    description = "G3-U70 wobble pair"
  ),
  det(
    "Ala",
    "Bacteria",
    4,
    "G",
    "acceptor_stem",
    "weak",
    pair_pos = 69L,
    pair_type = "wc",
    description = "G4-C69 base pair"
  ),
  det(
    "Ala",
    "Bacteria",
    69,
    "C",
    "acceptor_stem",
    "weak",
    pair_pos = 4L,
    pair_type = "wc",
    description = "G4-C69 base pair"
  ),
  det("Ala", "Bacteria", 20, "G", "d_arm", "weak", description = "G20"),
  det(
    "Ala",
    "Eukarya",
    73,
    "A",
    "discriminator",
    "strong",
    description = "A73 discriminator"
  ),
  det(
    "Ala",
    "Eukarya",
    3,
    "G",
    "acceptor_stem",
    "strong",
    pair_pos = 70L,
    pair_type = "wobble",
    universal = TRUE,
    description = "G3-U70 wobble pair"
  ),
  det(
    "Ala",
    "Eukarya",
    70,
    "U",
    "acceptor_stem",
    "strong",
    pair_pos = 3L,
    pair_type = "wobble",
    universal = TRUE,
    description = "G3-U70 wobble pair"
  ),

  # === Arg ===
  det(
    "Arg",
    "Bacteria",
    73,
    NA_character_,
    "discriminator",
    "strong",
    description = "A/G73 discriminator"
  ),
  det(
    "Arg",
    "Bacteria",
    2,
    "G",
    "acceptor_stem",
    "weak",
    pair_pos = 71L,
    pair_type = "wc",
    description = "G2-C71 base pair"
  ),
  det(
    "Arg",
    "Bacteria",
    71,
    "C",
    "acceptor_stem",
    "weak",
    pair_pos = 2L,
    pair_type = "wc",
    description = "G2-C71 base pair"
  ),
  det(
    "Arg",
    "Bacteria",
    3,
    "C",
    "acceptor_stem",
    "weak",
    pair_pos = 70L,
    pair_type = "wc",
    description = "C3-G70 base pair"
  ),
  det(
    "Arg",
    "Bacteria",
    70,
    "G",
    "acceptor_stem",
    "weak",
    pair_pos = 3L,
    pair_type = "wc",
    description = "C3-G70 base pair"
  ),
  det("Arg", "Bacteria", 20, "A", "d_arm", "weak", description = "A20"),
  det(
    "Arg",
    "Bacteria",
    35,
    "C",
    "anticodon_loop",
    "strong",
    universal = TRUE,
    description = "C35 anticodon"
  ),
  det(
    "Arg",
    "Bacteria",
    36,
    NA_character_,
    "anticodon_loop",
    "weak",
    description = "U/G36 anticodon"
  ),
  det(
    "Arg",
    "Eukarya",
    73,
    NA_character_,
    "discriminator",
    "strong",
    description = "A/C73 discriminator"
  ),
  det("Arg", "Eukarya", 20, "A", "d_arm", "weak", description = "A20"),
  det(
    "Arg",
    "Eukarya",
    35,
    "C",
    "anticodon_loop",
    "strong",
    universal = TRUE,
    description = "C35 anticodon"
  ),
  det(
    "Arg",
    "Eukarya",
    36,
    "G",
    "anticodon_loop",
    "weak",
    description = "G36 anticodon"
  ),

  # === Asn ===
  det(
    "Asn",
    "Bacteria",
    73,
    "G",
    "discriminator",
    "strong",
    description = "G73 discriminator"
  ),
  det(
    "Asn",
    "Bacteria",
    34,
    "G",
    "anticodon_loop",
    "strong",
    description = "G34 anticodon"
  ),
  det(
    "Asn",
    "Bacteria",
    1,
    NA_character_,
    "acceptor_stem",
    "weak",
    pair_pos = 72L,
    pair_type = "wc",
    description = "1-72 base pair"
  ),
  det(
    "Asn",
    "Bacteria",
    72,
    NA_character_,
    "acceptor_stem",
    "weak",
    pair_pos = 1L,
    pair_type = "wc",
    description = "1-72 base pair"
  ),
  det(
    "Asn",
    "Bacteria",
    35,
    "U",
    "anticodon_loop",
    "strong",
    universal = TRUE,
    description = "U35 anticodon"
  ),
  det(
    "Asn",
    "Bacteria",
    36,
    "U",
    "anticodon_loop",
    "strong",
    universal = TRUE,
    description = "U36 anticodon"
  ),
  det(
    "Asn",
    "Eukarya",
    73,
    "G",
    "discriminator",
    "strong",
    description = "G73 discriminator"
  ),
  det(
    "Asn",
    "Eukarya",
    35,
    "U",
    "anticodon_loop",
    "strong",
    universal = TRUE,
    description = "U35 anticodon"
  ),
  det(
    "Asn",
    "Eukarya",
    36,
    "U",
    "anticodon_loop",
    "strong",
    universal = TRUE,
    description = "U36 anticodon"
  ),

  # === Asp ===
  det(
    "Asp",
    "Bacteria",
    73,
    "G",
    "discriminator",
    "strong",
    universal = TRUE,
    description = "G73 discriminator"
  ),
  det("Asp", "Bacteria", 10, "G", "d_arm", "weak", description = "G10"),
  det(
    "Asp",
    "Bacteria",
    34,
    "G",
    "anticodon_loop",
    "strong",
    universal = TRUE,
    description = "G34 anticodon"
  ),
  det(
    "Asp",
    "Bacteria",
    35,
    "U",
    "anticodon_loop",
    "strong",
    universal = TRUE,
    description = "U35 anticodon"
  ),
  det(
    "Asp",
    "Bacteria",
    36,
    "C",
    "anticodon_loop",
    "strong",
    universal = TRUE,
    description = "C36 anticodon"
  ),
  det(
    "Asp",
    "Bacteria",
    38,
    "C",
    "anticodon_loop",
    "weak",
    description = "C38"
  ),
  det(
    "Asp",
    "Eukarya",
    73,
    "G",
    "discriminator",
    "strong",
    universal = TRUE,
    description = "G73 discriminator"
  ),
  det(
    "Asp",
    "Eukarya",
    34,
    "G",
    "anticodon_loop",
    "strong",
    universal = TRUE,
    description = "G34 anticodon"
  ),
  det(
    "Asp",
    "Eukarya",
    35,
    "U",
    "anticodon_loop",
    "strong",
    universal = TRUE,
    description = "U35 anticodon"
  ),
  det(
    "Asp",
    "Eukarya",
    36,
    "C",
    "anticodon_loop",
    "strong",
    universal = TRUE,
    description = "C36 anticodon"
  ),

  # === Cys ===
  det(
    "Cys",
    "Bacteria",
    73,
    "U",
    "discriminator",
    "strong",
    universal = TRUE,
    description = "U73 discriminator"
  ),
  det(
    "Cys",
    "Bacteria",
    2,
    "G",
    "acceptor_stem",
    "weak",
    pair_pos = 71L,
    pair_type = "wc",
    description = "G2-C71 base pair"
  ),
  det(
    "Cys",
    "Bacteria",
    71,
    "C",
    "acceptor_stem",
    "weak",
    pair_pos = 2L,
    pair_type = "wc",
    description = "G2-C71 base pair"
  ),
  det(
    "Cys",
    "Bacteria",
    3,
    "C",
    "acceptor_stem",
    "weak",
    pair_pos = 70L,
    pair_type = "wc",
    description = "C3-G70 base pair"
  ),
  det(
    "Cys",
    "Bacteria",
    70,
    "G",
    "acceptor_stem",
    "weak",
    pair_pos = 3L,
    pair_type = "wc",
    description = "C3-G70 base pair"
  ),
  det(
    "Cys",
    "Bacteria",
    34,
    "G",
    "anticodon_loop",
    "strong",
    universal = TRUE,
    description = "G34 anticodon"
  ),
  det(
    "Cys",
    "Bacteria",
    35,
    "C",
    "anticodon_loop",
    "strong",
    universal = TRUE,
    description = "C35 anticodon"
  ),
  det(
    "Cys",
    "Bacteria",
    36,
    "A",
    "anticodon_loop",
    "strong",
    universal = TRUE,
    description = "A36 anticodon"
  ),
  det(
    "Cys",
    "Eukarya",
    73,
    "U",
    "discriminator",
    "strong",
    universal = TRUE,
    description = "U73 discriminator"
  ),
  det(
    "Cys",
    "Eukarya",
    34,
    "G",
    "anticodon_loop",
    "strong",
    universal = TRUE,
    description = "G34 anticodon"
  ),
  det(
    "Cys",
    "Eukarya",
    35,
    "C",
    "anticodon_loop",
    "strong",
    universal = TRUE,
    description = "C35 anticodon"
  ),
  det(
    "Cys",
    "Eukarya",
    36,
    "A",
    "anticodon_loop",
    "strong",
    universal = TRUE,
    description = "A36 anticodon"
  ),

  # === Gln ===
  det(
    "Gln",
    "Bacteria",
    73,
    "G",
    "discriminator",
    "strong",
    description = "G73 discriminator"
  ),
  det(
    "Gln",
    "Bacteria",
    2,
    "G",
    "acceptor_stem",
    "strong",
    pair_pos = 71L,
    pair_type = "wc",
    description = "G2-C71 base pair"
  ),
  det(
    "Gln",
    "Bacteria",
    71,
    "C",
    "acceptor_stem",
    "strong",
    pair_pos = 2L,
    pair_type = "wc",
    description = "G2-C71 base pair"
  ),
  det(
    "Gln",
    "Bacteria",
    3,
    "G",
    "acceptor_stem",
    "strong",
    pair_pos = 70L,
    pair_type = "wc",
    description = "G3-C70 base pair"
  ),
  det(
    "Gln",
    "Bacteria",
    70,
    "C",
    "acceptor_stem",
    "strong",
    pair_pos = 3L,
    pair_type = "wc",
    description = "G3-C70 base pair"
  ),
  det("Gln", "Bacteria", 10, "U", "d_arm", "weak", description = "U10"),
  det(
    "Gln",
    "Bacteria",
    34,
    NA_character_,
    "anticodon_loop",
    "strong",
    universal = TRUE,
    description = "anticodon pos 34"
  ),
  det(
    "Gln",
    "Bacteria",
    35,
    "U",
    "anticodon_loop",
    "strong",
    universal = TRUE,
    description = "U35 anticodon"
  ),
  det(
    "Gln",
    "Bacteria",
    36,
    "G",
    "anticodon_loop",
    "strong",
    universal = TRUE,
    description = "G36 anticodon"
  ),
  det(
    "Gln",
    "Eukarya",
    73,
    NA_character_,
    "discriminator",
    "strong",
    description = "discriminator base 73"
  ),
  det(
    "Gln",
    "Eukarya",
    2,
    NA_character_,
    "acceptor_stem",
    "weak",
    pair_pos = 71L,
    pair_type = "wc",
    description = "2-71 base pair"
  ),
  det(
    "Gln",
    "Eukarya",
    71,
    NA_character_,
    "acceptor_stem",
    "weak",
    pair_pos = 2L,
    pair_type = "wc",
    description = "2-71 base pair"
  ),
  det(
    "Gln",
    "Eukarya",
    34,
    NA_character_,
    "anticodon_loop",
    "strong",
    universal = TRUE,
    description = "anticodon pos 34"
  ),
  det(
    "Gln",
    "Eukarya",
    35,
    "U",
    "anticodon_loop",
    "strong",
    universal = TRUE,
    description = "U35 anticodon"
  ),
  det(
    "Gln",
    "Eukarya",
    36,
    "G",
    "anticodon_loop",
    "strong",
    universal = TRUE,
    description = "G36 anticodon"
  ),

  # === Glu ===
  # Note: position 73 is NOT a GluRS identity determinant in E. coli
  # (Glu and Thr are the two stated exceptions in Giege & Eriani 2023).
  det(
    "Glu",
    "Bacteria",
    1,
    NA_character_,
    "acceptor_stem",
    "weak",
    pair_pos = 72L,
    pair_type = "wc",
    description = "1-72 base pair"
  ),
  det(
    "Glu",
    "Bacteria",
    72,
    NA_character_,
    "acceptor_stem",
    "weak",
    pair_pos = 1L,
    pair_type = "wc",
    description = "1-72 base pair"
  ),
  det(
    "Glu",
    "Bacteria",
    34,
    NA_character_,
    "anticodon_loop",
    "strong",
    universal = TRUE,
    description = "anticodon pos 34"
  ),
  det(
    "Glu",
    "Bacteria",
    35,
    "U",
    "anticodon_loop",
    "strong",
    universal = TRUE,
    description = "U35 anticodon"
  ),
  det(
    "Glu",
    "Bacteria",
    36,
    "C",
    "anticodon_loop",
    "strong",
    universal = TRUE,
    description = "C36 anticodon"
  ),
  det(
    "Glu",
    "Eukarya",
    73,
    NA_character_,
    "discriminator",
    "strong",
    description = "discriminator base 73"
  ),
  det(
    "Glu",
    "Eukarya",
    34,
    NA_character_,
    "anticodon_loop",
    "strong",
    universal = TRUE,
    description = "anticodon pos 34"
  ),
  det(
    "Glu",
    "Eukarya",
    35,
    "U",
    "anticodon_loop",
    "strong",
    universal = TRUE,
    description = "U35 anticodon"
  ),
  det(
    "Glu",
    "Eukarya",
    36,
    "C",
    "anticodon_loop",
    "strong",
    universal = TRUE,
    description = "C36 anticodon"
  ),

  # === Gly ===
  det(
    "Gly",
    "Bacteria",
    73,
    "U",
    "discriminator",
    "strong",
    description = "U73 discriminator"
  ),
  det(
    "Gly",
    "Bacteria",
    1,
    "G",
    "acceptor_stem",
    "weak",
    pair_pos = 72L,
    pair_type = "wc",
    description = "G1-C72 base pair"
  ),
  det(
    "Gly",
    "Bacteria",
    72,
    "C",
    "acceptor_stem",
    "weak",
    pair_pos = 1L,
    pair_type = "wc",
    description = "G1-C72 base pair"
  ),
  det(
    "Gly",
    "Bacteria",
    2,
    "C",
    "acceptor_stem",
    "weak",
    pair_pos = 71L,
    pair_type = "wc",
    description = "C2-G71 base pair"
  ),
  det(
    "Gly",
    "Bacteria",
    71,
    "G",
    "acceptor_stem",
    "weak",
    pair_pos = 2L,
    pair_type = "wc",
    description = "C2-G71 base pair"
  ),
  det(
    "Gly",
    "Bacteria",
    3,
    "G",
    "acceptor_stem",
    "weak",
    pair_pos = 70L,
    pair_type = "wc",
    description = "G3-C70 base pair"
  ),
  det(
    "Gly",
    "Bacteria",
    70,
    "C",
    "acceptor_stem",
    "weak",
    pair_pos = 3L,
    pair_type = "wc",
    description = "G3-C70 base pair"
  ),
  det("Gly", "Bacteria", 10, "G", "d_arm", "weak", description = "G10"),
  det(
    "Gly",
    "Bacteria",
    34,
    NA_character_,
    "anticodon_loop",
    "strong",
    universal = TRUE,
    description = "anticodon pos 34"
  ),
  det(
    "Gly",
    "Bacteria",
    35,
    "C",
    "anticodon_loop",
    "strong",
    universal = TRUE,
    description = "C35 anticodon"
  ),
  det(
    "Gly",
    "Bacteria",
    36,
    "C",
    "anticodon_loop",
    "strong",
    universal = TRUE,
    description = "C36 anticodon"
  ),
  det(
    "Gly",
    "Eukarya",
    73,
    "A",
    "discriminator",
    "strong",
    description = "A73 discriminator"
  ),
  det(
    "Gly",
    "Eukarya",
    34,
    NA_character_,
    "anticodon_loop",
    "strong",
    universal = TRUE,
    description = "anticodon pos 34"
  ),
  det(
    "Gly",
    "Eukarya",
    35,
    "C",
    "anticodon_loop",
    "strong",
    universal = TRUE,
    description = "C35 anticodon"
  ),
  det(
    "Gly",
    "Eukarya",
    36,
    "C",
    "anticodon_loop",
    "strong",
    universal = TRUE,
    description = "C36 anticodon"
  ),

  # === His ===
  det(
    "His",
    "Bacteria",
    -1,
    "G",
    "acceptor_stem",
    "strong",
    universal = TRUE,
    description = "G-1 extra 5' nucleotide"
  ),
  det(
    "His",
    "Bacteria",
    73,
    "C",
    "discriminator",
    "strong",
    universal = TRUE,
    description = "C73 discriminator"
  ),
  det(
    "His",
    "Bacteria",
    34,
    "G",
    "anticodon_loop",
    "strong",
    description = "G34 anticodon"
  ),
  det(
    "His",
    "Bacteria",
    35,
    "U",
    "anticodon_loop",
    "strong",
    description = "U35 anticodon"
  ),
  det(
    "His",
    "Bacteria",
    36,
    "G",
    "anticodon_loop",
    "strong",
    description = "G36 anticodon"
  ),
  det(
    "His",
    "Eukarya",
    -1,
    "G",
    "acceptor_stem",
    "strong",
    universal = TRUE,
    description = "G-1 extra 5' nucleotide"
  ),
  det(
    "His",
    "Eukarya",
    73,
    "C",
    "discriminator",
    "strong",
    universal = TRUE,
    description = "C73 discriminator"
  ),
  det(
    "His",
    "Eukarya",
    34,
    "G",
    "anticodon_loop",
    "weak",
    description = "G34 anticodon"
  ),
  det(
    "His",
    "Eukarya",
    35,
    "U",
    "anticodon_loop",
    "weak",
    description = "U35 anticodon"
  ),
  det(
    "His",
    "Eukarya",
    36,
    "G",
    "anticodon_loop",
    "weak",
    description = "G36 anticodon"
  ),

  # === Ile ===
  det(
    "Ile",
    "Bacteria",
    73,
    "A",
    "discriminator",
    "strong",
    description = "A73 discriminator"
  ),
  det(
    "Ile",
    "Bacteria",
    34,
    NA_character_,
    "anticodon_loop",
    "strong",
    universal = TRUE,
    description = "anticodon pos 34 (k2C34 modification)"
  ),
  det(
    "Ile",
    "Bacteria",
    35,
    "A",
    "anticodon_loop",
    "strong",
    universal = TRUE,
    description = "A35 anticodon"
  ),
  det(
    "Ile",
    "Bacteria",
    36,
    "U",
    "anticodon_loop",
    "strong",
    universal = TRUE,
    description = "U36 anticodon"
  ),
  det(
    "Ile",
    "Bacteria",
    12,
    NA_character_,
    "d_arm",
    "weak",
    description = "position 12"
  ),
  det(
    "Ile",
    "Eukarya",
    73,
    "A",
    "discriminator",
    "strong",
    description = "A73 discriminator"
  ),
  det(
    "Ile",
    "Eukarya",
    34,
    NA_character_,
    "anticodon_loop",
    "strong",
    universal = TRUE,
    description = "anticodon pos 34"
  ),
  det(
    "Ile",
    "Eukarya",
    35,
    "A",
    "anticodon_loop",
    "strong",
    universal = TRUE,
    description = "A35 anticodon"
  ),
  det(
    "Ile",
    "Eukarya",
    36,
    "U",
    "anticodon_loop",
    "strong",
    universal = TRUE,
    description = "U36 anticodon"
  ),

  # === Leu ===
  det(
    "Leu",
    "Bacteria",
    73,
    "A",
    "discriminator",
    "strong",
    description = "A73 discriminator"
  ),
  det(
    "Leu",
    "Bacteria",
    35,
    "A",
    "anticodon_loop",
    "strong",
    description = "A35 anticodon"
  ),
  det(
    "Leu",
    "Bacteria",
    NA_integer_,
    NA_character_,
    "variable_loop",
    "strong",
    description = "long variable arm (structural determinant)"
  ),
  det(
    "Leu",
    "Eukarya",
    73,
    "A",
    "discriminator",
    "strong",
    description = "A73 discriminator"
  ),
  det(
    "Leu",
    "Eukarya",
    35,
    NA_character_,
    "anticodon_loop",
    "strong",
    description = "anticodon pos 35"
  ),
  det(
    "Leu",
    "Eukarya",
    NA_integer_,
    NA_character_,
    "variable_loop",
    "strong",
    description = "long variable arm (structural determinant)"
  ),

  # === Lys ===
  det(
    "Lys",
    "Bacteria",
    73,
    "A",
    "discriminator",
    "strong",
    description = "A73 discriminator"
  ),
  det(
    "Lys",
    "Bacteria",
    34,
    NA_character_,
    "anticodon_loop",
    "strong",
    universal = TRUE,
    description = "anticodon pos 34"
  ),
  det(
    "Lys",
    "Bacteria",
    35,
    "U",
    "anticodon_loop",
    "strong",
    universal = TRUE,
    description = "U35 anticodon"
  ),
  det(
    "Lys",
    "Bacteria",
    36,
    "U",
    "anticodon_loop",
    "strong",
    universal = TRUE,
    description = "U36 anticodon"
  ),
  det(
    "Lys",
    "Eukarya",
    73,
    NA_character_,
    "discriminator",
    "strong",
    description = "discriminator base 73"
  ),
  det(
    "Lys",
    "Eukarya",
    34,
    NA_character_,
    "anticodon_loop",
    "strong",
    universal = TRUE,
    description = "anticodon pos 34"
  ),
  det(
    "Lys",
    "Eukarya",
    35,
    "U",
    "anticodon_loop",
    "strong",
    universal = TRUE,
    description = "U35 anticodon"
  ),
  det(
    "Lys",
    "Eukarya",
    36,
    "U",
    "anticodon_loop",
    "strong",
    universal = TRUE,
    description = "U36 anticodon"
  ),

  # === Met ===
  det(
    "Met",
    "Bacteria",
    73,
    "A",
    "discriminator",
    "strong",
    description = "A73 discriminator"
  ),
  det(
    "Met",
    "Bacteria",
    1,
    "G",
    "acceptor_stem",
    "weak",
    pair_pos = 72L,
    pair_type = "wc",
    description = "G1-C72 base pair"
  ),
  det(
    "Met",
    "Bacteria",
    72,
    "C",
    "acceptor_stem",
    "weak",
    pair_pos = 1L,
    pair_type = "wc",
    description = "G1-C72 base pair"
  ),
  det(
    "Met",
    "Bacteria",
    2,
    "C",
    "acceptor_stem",
    "weak",
    pair_pos = 71L,
    pair_type = "wc",
    description = "C2-G71 base pair"
  ),
  det(
    "Met",
    "Bacteria",
    71,
    "G",
    "acceptor_stem",
    "weak",
    pair_pos = 2L,
    pair_type = "wc",
    description = "C2-G71 base pair"
  ),
  det(
    "Met",
    "Bacteria",
    34,
    "C",
    "anticodon_loop",
    "strong",
    universal = TRUE,
    description = "C34 anticodon"
  ),
  det(
    "Met",
    "Bacteria",
    35,
    "A",
    "anticodon_loop",
    "strong",
    universal = TRUE,
    description = "A35 anticodon"
  ),
  det(
    "Met",
    "Bacteria",
    36,
    "U",
    "anticodon_loop",
    "strong",
    universal = TRUE,
    description = "U36 anticodon"
  ),
  det(
    "Met",
    "Eukarya",
    73,
    "A",
    "discriminator",
    "strong",
    description = "A73 discriminator"
  ),
  det(
    "Met",
    "Eukarya",
    34,
    "C",
    "anticodon_loop",
    "strong",
    universal = TRUE,
    description = "C34 anticodon"
  ),
  det(
    "Met",
    "Eukarya",
    35,
    "A",
    "anticodon_loop",
    "strong",
    universal = TRUE,
    description = "A35 anticodon"
  ),
  det(
    "Met",
    "Eukarya",
    36,
    "U",
    "anticodon_loop",
    "strong",
    universal = TRUE,
    description = "U36 anticodon"
  ),

  # === Phe ===
  det(
    "Phe",
    "Bacteria",
    73,
    "A",
    "discriminator",
    "strong",
    description = "A73 discriminator"
  ),
  det(
    "Phe",
    "Bacteria",
    10,
    "C",
    "d_arm",
    "weak",
    pair_pos = 25L,
    pair_type = "wc",
    description = "C10-G25 base pair"
  ),
  det(
    "Phe",
    "Bacteria",
    25,
    "G",
    "d_arm",
    "weak",
    pair_pos = 10L,
    pair_type = "wc",
    description = "C10-G25 base pair"
  ),
  det("Phe", "Bacteria", 20, "U", "d_arm", "weak", description = "U20"),
  det("Phe", "Bacteria", 45, "U", "variable_loop", "weak", description = "U45"),
  det("Phe", "Bacteria", 59, "U", "t_arm", "weak", description = "U59"),
  det(
    "Phe",
    "Bacteria",
    34,
    "G",
    "anticodon_loop",
    "strong",
    universal = TRUE,
    description = "G34 anticodon"
  ),
  det(
    "Phe",
    "Bacteria",
    35,
    "A",
    "anticodon_loop",
    "strong",
    universal = TRUE,
    description = "A35 anticodon"
  ),
  det(
    "Phe",
    "Bacteria",
    36,
    "A",
    "anticodon_loop",
    "strong",
    universal = TRUE,
    description = "A36 anticodon"
  ),
  det(
    "Phe",
    "Eukarya",
    73,
    "A",
    "discriminator",
    "strong",
    description = "A73 discriminator"
  ),
  det(
    "Phe",
    "Eukarya",
    20,
    NA_character_,
    "d_arm",
    "weak",
    description = "position 20"
  ),
  det(
    "Phe",
    "Eukarya",
    34,
    "G",
    "anticodon_loop",
    "strong",
    universal = TRUE,
    description = "G34 anticodon"
  ),
  det(
    "Phe",
    "Eukarya",
    35,
    "A",
    "anticodon_loop",
    "strong",
    universal = TRUE,
    description = "A35 anticodon"
  ),
  det(
    "Phe",
    "Eukarya",
    36,
    "A",
    "anticodon_loop",
    "strong",
    universal = TRUE,
    description = "A36 anticodon"
  ),

  # === Pro ===
  det(
    "Pro",
    "Bacteria",
    73,
    "A",
    "discriminator",
    "strong",
    description = "A73 discriminator"
  ),
  det("Pro", "Bacteria", 72, "G", "acceptor_stem", "weak", description = "G72"),
  det(
    "Pro",
    "Bacteria",
    34,
    NA_character_,
    "anticodon_loop",
    "strong",
    universal = TRUE,
    description = "anticodon pos 34"
  ),
  det(
    "Pro",
    "Bacteria",
    35,
    "G",
    "anticodon_loop",
    "strong",
    universal = TRUE,
    description = "G35 anticodon"
  ),
  det(
    "Pro",
    "Bacteria",
    36,
    "G",
    "anticodon_loop",
    "strong",
    universal = TRUE,
    description = "G36 anticodon"
  ),
  det(
    "Pro",
    "Eukarya",
    73,
    "A",
    "discriminator",
    "strong",
    description = "A73 discriminator"
  ),
  det(
    "Pro",
    "Eukarya",
    34,
    NA_character_,
    "anticodon_loop",
    "strong",
    universal = TRUE,
    description = "anticodon pos 34"
  ),
  det(
    "Pro",
    "Eukarya",
    35,
    "G",
    "anticodon_loop",
    "strong",
    universal = TRUE,
    description = "G35 anticodon"
  ),
  det(
    "Pro",
    "Eukarya",
    36,
    "G",
    "anticodon_loop",
    "strong",
    universal = TRUE,
    description = "G36 anticodon"
  ),

  # === Ser ===
  det(
    "Ser",
    "Bacteria",
    73,
    "G",
    "discriminator",
    "strong",
    description = "G73 discriminator"
  ),
  det(
    "Ser",
    "Bacteria",
    NA_integer_,
    NA_character_,
    "variable_loop",
    "strong",
    description = "long variable arm (structural determinant)"
  ),
  det(
    "Ser",
    "Bacteria",
    12,
    "G",
    "d_arm",
    "weak",
    pair_pos = 23L,
    pair_type = "wc",
    description = "G12-C23 base pair"
  ),
  det(
    "Ser",
    "Bacteria",
    23,
    "C",
    "d_arm",
    "weak",
    pair_pos = 12L,
    pair_type = "wc",
    description = "G12-C23 base pair"
  ),
  det(
    "Ser",
    "Eukarya",
    73,
    "G",
    "discriminator",
    "strong",
    description = "G73 discriminator"
  ),
  det(
    "Ser",
    "Eukarya",
    NA_integer_,
    NA_character_,
    "variable_loop",
    "strong",
    description = "long variable arm (structural determinant)"
  ),

  # === Thr ===
  # Note: position 73 is NOT a ThrRS identity determinant in E. coli
  # (Glu and Thr are the two stated exceptions in Giege & Eriani 2023;
  # the E. coli tRNA-Thr discriminator is A73, shared with Val).
  det(
    "Thr",
    "Bacteria",
    2,
    "G",
    "acceptor_stem",
    "strong",
    pair_pos = 71L,
    pair_type = "wc",
    description = "G2-C71 base pair"
  ),
  det(
    "Thr",
    "Bacteria",
    71,
    "C",
    "acceptor_stem",
    "strong",
    pair_pos = 2L,
    pair_type = "wc",
    description = "G2-C71 base pair"
  ),
  det(
    "Thr",
    "Bacteria",
    3,
    NA_character_,
    "acceptor_stem",
    "weak",
    pair_pos = 70L,
    pair_type = "wc",
    description = "3-70 base pair"
  ),
  det(
    "Thr",
    "Bacteria",
    70,
    NA_character_,
    "acceptor_stem",
    "weak",
    pair_pos = 3L,
    pair_type = "wc",
    description = "3-70 base pair"
  ),
  det(
    "Thr",
    "Bacteria",
    35,
    "G",
    "anticodon_loop",
    "strong",
    universal = TRUE,
    description = "G35 anticodon"
  ),
  det(
    "Thr",
    "Bacteria",
    36,
    "U",
    "anticodon_loop",
    "strong",
    universal = TRUE,
    description = "U36 anticodon"
  ),
  det(
    "Thr",
    "Eukarya",
    73,
    NA_character_,
    "discriminator",
    "strong",
    description = "discriminator base 73"
  ),
  det(
    "Thr",
    "Eukarya",
    2,
    NA_character_,
    "acceptor_stem",
    "weak",
    pair_pos = 71L,
    pair_type = "wc",
    description = "2-71 base pair"
  ),
  det(
    "Thr",
    "Eukarya",
    71,
    NA_character_,
    "acceptor_stem",
    "weak",
    pair_pos = 2L,
    pair_type = "wc",
    description = "2-71 base pair"
  ),
  det(
    "Thr",
    "Eukarya",
    35,
    "G",
    "anticodon_loop",
    "strong",
    universal = TRUE,
    description = "G35 anticodon"
  ),
  det(
    "Thr",
    "Eukarya",
    36,
    "U",
    "anticodon_loop",
    "strong",
    universal = TRUE,
    description = "U36 anticodon"
  ),

  # === Trp ===
  det(
    "Trp",
    "Bacteria",
    73,
    "A",
    "discriminator",
    "strong",
    description = "A73 discriminator"
  ),
  det(
    "Trp",
    "Bacteria",
    1,
    "G",
    "acceptor_stem",
    "weak",
    pair_pos = 72L,
    pair_type = "wc",
    description = "G1-C72 base pair"
  ),
  det(
    "Trp",
    "Bacteria",
    72,
    "C",
    "acceptor_stem",
    "weak",
    pair_pos = 1L,
    pair_type = "wc",
    description = "G1-C72 base pair"
  ),
  det(
    "Trp",
    "Bacteria",
    5,
    "G",
    "acceptor_stem",
    "weak",
    pair_pos = 68L,
    pair_type = "wc",
    description = "G5-C68 base pair"
  ),
  det(
    "Trp",
    "Bacteria",
    68,
    "C",
    "acceptor_stem",
    "weak",
    pair_pos = 5L,
    pair_type = "wc",
    description = "G5-C68 base pair"
  ),
  det(
    "Trp",
    "Bacteria",
    34,
    "C",
    "anticodon_loop",
    "strong",
    universal = TRUE,
    description = "C34 anticodon"
  ),
  det(
    "Trp",
    "Bacteria",
    35,
    "C",
    "anticodon_loop",
    "strong",
    universal = TRUE,
    description = "C35 anticodon"
  ),
  det(
    "Trp",
    "Bacteria",
    36,
    "A",
    "anticodon_loop",
    "strong",
    universal = TRUE,
    description = "A36 anticodon"
  ),
  det(
    "Trp",
    "Eukarya",
    73,
    "G",
    "discriminator",
    "strong",
    description = "G73 discriminator"
  ),
  det(
    "Trp",
    "Eukarya",
    1,
    "C",
    "acceptor_stem",
    "weak",
    pair_pos = 72L,
    pair_type = "wc",
    description = "C1-G72 base pair"
  ),
  det(
    "Trp",
    "Eukarya",
    72,
    "G",
    "acceptor_stem",
    "weak",
    pair_pos = 1L,
    pair_type = "wc",
    description = "C1-G72 base pair"
  ),
  det(
    "Trp",
    "Eukarya",
    34,
    "C",
    "anticodon_loop",
    "strong",
    universal = TRUE,
    description = "C34 anticodon"
  ),
  det(
    "Trp",
    "Eukarya",
    35,
    "C",
    "anticodon_loop",
    "strong",
    universal = TRUE,
    description = "C35 anticodon"
  ),
  det(
    "Trp",
    "Eukarya",
    36,
    "A",
    "anticodon_loop",
    "strong",
    universal = TRUE,
    description = "A36 anticodon"
  ),

  # === Tyr ===
  det(
    "Tyr",
    "Bacteria",
    73,
    "A",
    "discriminator",
    "strong",
    description = "A73 discriminator"
  ),
  det(
    "Tyr",
    "Bacteria",
    1,
    "G",
    "acceptor_stem",
    "strong",
    pair_pos = 72L,
    pair_type = "wc",
    description = "G1-C72 base pair"
  ),
  det(
    "Tyr",
    "Bacteria",
    72,
    "C",
    "acceptor_stem",
    "strong",
    pair_pos = 1L,
    pair_type = "wc",
    description = "G1-C72 base pair"
  ),
  det(
    "Tyr",
    "Bacteria",
    34,
    "G",
    "anticodon_loop",
    "strong",
    universal = TRUE,
    description = "G34 anticodon"
  ),
  det(
    "Tyr",
    "Bacteria",
    35,
    "U",
    "anticodon_loop",
    "strong",
    universal = TRUE,
    description = "U35 anticodon"
  ),
  det(
    "Tyr",
    "Bacteria",
    36,
    "A",
    "anticodon_loop",
    "strong",
    universal = TRUE,
    description = "A36 anticodon"
  ),
  det(
    "Tyr",
    "Eukarya",
    73,
    "A",
    "discriminator",
    "strong",
    description = "A73 discriminator"
  ),
  det(
    "Tyr",
    "Eukarya",
    1,
    "G",
    "acceptor_stem",
    "strong",
    pair_pos = 72L,
    pair_type = "wc",
    description = "G1-C72 base pair"
  ),
  det(
    "Tyr",
    "Eukarya",
    72,
    "C",
    "acceptor_stem",
    "strong",
    pair_pos = 1L,
    pair_type = "wc",
    description = "G1-C72 base pair"
  ),
  det(
    "Tyr",
    "Eukarya",
    34,
    "G",
    "anticodon_loop",
    "strong",
    universal = TRUE,
    description = "G34 anticodon"
  ),
  det(
    "Tyr",
    "Eukarya",
    35,
    "U",
    "anticodon_loop",
    "strong",
    universal = TRUE,
    description = "U35 anticodon"
  ),
  det(
    "Tyr",
    "Eukarya",
    36,
    "A",
    "anticodon_loop",
    "strong",
    universal = TRUE,
    description = "A36 anticodon"
  ),

  # === Val ===
  det(
    "Val",
    "Bacteria",
    73,
    "A",
    "discriminator",
    "strong",
    description = "A73 discriminator"
  ),
  det(
    "Val",
    "Bacteria",
    34,
    NA_character_,
    "anticodon_loop",
    "strong",
    universal = TRUE,
    description = "anticodon pos 34"
  ),
  det(
    "Val",
    "Bacteria",
    35,
    "A",
    "anticodon_loop",
    "strong",
    universal = TRUE,
    description = "A35 anticodon"
  ),
  det(
    "Val",
    "Bacteria",
    36,
    "C",
    "anticodon_loop",
    "strong",
    universal = TRUE,
    description = "C36 anticodon"
  ),
  det("Val", "Bacteria", 20, "A", "d_arm", "weak", description = "A20"),
  det(
    "Val",
    "Eukarya",
    73,
    "A",
    "discriminator",
    "strong",
    description = "A73 discriminator"
  ),
  det(
    "Val",
    "Eukarya",
    34,
    NA_character_,
    "anticodon_loop",
    "strong",
    universal = TRUE,
    description = "anticodon pos 34"
  ),
  det(
    "Val",
    "Eukarya",
    35,
    "A",
    "anticodon_loop",
    "strong",
    universal = TRUE,
    description = "A35 anticodon"
  ),
  det(
    "Val",
    "Eukarya",
    36,
    "C",
    "anticodon_loop",
    "strong",
    universal = TRUE,
    description = "C36 anticodon"
  )
)

# -- Antideterminants (Table 3) --
# Elements that prevent aminoacylation by a specific aaRS

antideterminants <- bind_rows(
  # E. coli antideterminants
  anti(
    "Trp",
    "Bacteria",
    1,
    "A",
    "acceptor_stem",
    against_aars = "MetRS",
    pair_pos = 72L,
    pair_type = "wc",
    description = "A1-U72 in tRNA-Trp blocks MetRS"
  ),
  anti(
    "Trp",
    "Bacteria",
    72,
    "U",
    "acceptor_stem",
    against_aars = "MetRS",
    pair_pos = 1L,
    pair_type = "wc",
    description = "A1-U72 in tRNA-Trp blocks MetRS"
  ),
  anti(
    "Leu",
    "Bacteria",
    2,
    "C",
    "acceptor_stem",
    against_aars = "SerRS",
    pair_pos = 71L,
    pair_type = "wc",
    description = "C2-G71 in tRNA-Leu blocks SerRS"
  ),
  anti(
    "Leu",
    "Bacteria",
    71,
    "G",
    "acceptor_stem",
    against_aars = "SerRS",
    pair_pos = 2L,
    pair_type = "wc",
    description = "C2-G71 in tRNA-Leu blocks SerRS"
  ),
  anti(
    "Ile",
    "Bacteria",
    34,
    "C",
    "anticodon_loop",
    against_aars = "MetRS",
    description = "k2C34 in tRNA-Ile blocks MetRS (lysidine modification)"
  ),
  anti(
    "Ser",
    "Bacteria",
    73,
    "G",
    "discriminator",
    against_aars = "LeuRS",
    description = "G73 in tRNA-Ser blocks LeuRS"
  ),
  anti(
    "Ser",
    "Bacteria",
    73,
    "G",
    "discriminator",
    against_aars = "TyrRS",
    description = "G73 in tRNA-Ser blocks TyrRS"
  ),

  # S. cerevisiae / H. sapiens antideterminants
  anti(
    "Ala",
    "Eukarya",
    3,
    "G",
    "acceptor_stem",
    against_aars = "ThrRS",
    pair_pos = 70L,
    pair_type = "wobble",
    description = "G3-U70 in tRNA-Ala blocks ThrRS"
  ),
  anti(
    "Ala",
    "Eukarya",
    70,
    "U",
    "acceptor_stem",
    against_aars = "ThrRS",
    pair_pos = 3L,
    pair_type = "wobble",
    description = "G3-U70 in tRNA-Ala blocks ThrRS"
  ),
  anti(
    "Leu",
    "Eukarya",
    73,
    "A",
    "discriminator",
    against_aars = "SerRS",
    description = "A73 in tRNA-Leu blocks SerRS"
  )
)

# -- Organism map --
organism_map <- tibble(
  organism = c(
    "Escherichia coli",
    "Saccharomyces cerevisiae",
    "Homo sapiens"
  ),
  domain = c("Bacteria", "Eukarya", "Eukarya")
)

# -- Save RDS files --
cli::cli_inform("Saving determinants ({nrow(determinants)} rows).")
saveRDS(determinants, file.path(out_dir, "determinants.rds"))

cli::cli_inform("Saving antideterminants ({nrow(antideterminants)} rows).")
saveRDS(antideterminants, file.path(out_dir, "antideterminants.rds"))

cli::cli_inform("Saving organism map ({nrow(organism_map)} rows).")
saveRDS(organism_map, file.path(out_dir, "organism_map.rds"))

cli::cli_inform("Done.")
