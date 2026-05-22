# match_modomics_to_refs warns when no alignment passes min_identity

    Code
      result <- match_modomics_to_refs(modomics_entries, ref_seq, min_identity = 1.01)
    Condition
      Warning:
      ! No reference matched MODOMICS "Ala-AGC" (length 13) at min_identity = 1.01; best identity was 1.

