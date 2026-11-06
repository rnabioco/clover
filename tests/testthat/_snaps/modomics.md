# match_modomics_to_refs warns when no alignment passes min_identity

    Code
      result <- match_modomics_to_refs(modomics_entries, ref_seq, min_identity = 1.01)
    Condition
      Warning:
      ! No reference matched MODOMICS "Ala-AGC" (length 13) at min_identity = 1.01 and min_coverage = 0.8; best identity was 1 and best coverage was 0.08.

# match_modomics_to_refs rejects short local fragments

    Code
      result <- match_modomics_to_refs(list(entry), ref_seq, min_identity = 0.7)
    Condition
      Warning:
      ! No reference matched MODOMICS "Ala-AGC" (length 40) at min_identity = 0.7 and min_coverage = 0.8; best identity was 1 and best coverage was 0.25.

