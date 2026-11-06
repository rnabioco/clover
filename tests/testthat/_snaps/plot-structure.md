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

# convert_sprinzl_positions warns on unmatched and drops rows

    Code
      result <- convert_sprinzl_positions(df, "pos", trna_coords)
    Condition
      Warning:
      Sprinzl position "99" not found; dropping 1 row.

# plot_tRNA_structure errors when tRNA not in sprinzl_coords

    Code
      plot_tRNA_structure(trna, org, sprinzl_coords = fake_coords)
    Condition
      Error in `plot_tRNA_structure()`:
      ! Could not find "tRNA-Ala-GGC" in `sprinzl_coords`.

