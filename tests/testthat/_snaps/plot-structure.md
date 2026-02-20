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

