# list_pipeline_files validates types

    Code
      list_pipeline_files(config, types = "invalid_type")
    Condition
      Error in `match.arg()`:
      ! 'arg' should be one of "charging", "bcerror", "odds_ratios", "align_stats"

