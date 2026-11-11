# structure_trnas errors for unknown organism

    Code
      structure_trnas("Nonexistent organism")
    Condition
      Error in `structure_org_dir()`:
      ! No structure data found for "Nonexistent organism".
      i Use `structure_organisms()` to list available organisms.

# plot_tRNA_structure errors for unknown organism

    Code
      plot_tRNA_structure("tRNA-Ala-GGC", "Nonexistent organism")
    Condition
      Error in `structure_org_dir()`:
      ! No structure data found for "Nonexistent organism".
      i Use `structure_organisms()` to list available organisms.

# plot_tRNA_structure errors for unknown tRNA

    Code
      plot_tRNA_structure("tRNA-Fake-XXX", org)
    Condition
      Error in `plot_tRNA_structure()`:
      ! No structure SVG found for "tRNA-Fake-XXX".
      i Use `structure_trnas()` to list available tRNAs.

# structure_to_png errors for missing file

    Code
      structure_to_png("nonexistent.svg")
    Condition
      Error in `structure_to_png()`:
      ! SVG file not found: 'nonexistent.svg'.

# structure_html errors on missing file

    Code
      structure_html("nonexistent.svg")
    Condition
      Error in `structure_html()`:
      ! SVG file not found: 'nonexistent.svg'.

# plot_tRNA_diff errors when the reference does not match

    Code
      plot_tRNA_diff("tRNA-Gly-TCC", "Escherichia coli", reference_coords = ecoli,
        reference_id = "tRNA-Gly-GCC-1-1", query_coords = ecoli, query_id = "tRNA-Gly-UCC-1-1")
    Condition
      Error in `check_diff_matches_structure()`:
      ! `reference_coords` does not match the "tRNA-Gly-TCC" structure.
      x 43 of 75 positions differ between "tRNA-Gly-GCC-1-1" and the bundled cloverleaf.
      i Check that `reference_id` is the same tRNA as `trna`.

