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
      ! No cloverleaf structure SVG found for "tRNA-Fake-XXX".
      i Use `structure_trnas("Escherichia coli", layout = "cloverleaf")` to list available tRNAs.

# structure_html errors on missing file

    Code
      structure_html("nonexistent.svg")
    Condition
      Error in `structure_html()`:
      ! SVG file not found: 'nonexistent.svg'.

# plot_tRNA_structure(layout = 'elbow') errors clearly when missing

    Code
      plot_tRNA_structure("tRNA-Leu-CAA", "Escherichia coli", layout = "elbow")
    Condition
      Error in `plot_tRNA_structure()`:
      ! No elbow structure SVG found for "tRNA-Leu-CAA".
      i Use `structure_trnas("Escherichia coli", layout = "elbow")` to list available tRNAs.

