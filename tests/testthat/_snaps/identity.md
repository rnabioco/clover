# identity_elements errors for unsupported organism

    Code
      identity_elements("Mus musculus")
    Condition
      Error in `identity_elements()`:
      ! Organism "Mus musculus" is not supported.
      i Supported organisms: 'Escherichia coli', 'Saccharomyces cerevisiae', 'Homo sapiens'.
      i Use `identity_organisms()` to list them.

# map_identity_to_trna errors for missing tRNA

    Code
      map_identity_to_trna(elements, mock_coords, "tRNA-Ala-GGC")
    Condition
      Error in `map_identity_to_trna()`:
      ! tRNA "tRNA-Ala-GGC" not found in `sprinzl_coords`.

