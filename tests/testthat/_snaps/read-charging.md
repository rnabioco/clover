# read_charging_multi errors on unnamed paths

    Code
      read_charging_multi(c("/fake/path1", "/fake/path2"))
    Condition
      Error in `read_multi()`:
      ! `paths` must be a named character vector.

# read_odds_ratios_multi errors on unnamed paths

    Code
      read_odds_ratios_multi(c("/fake/path1"))
    Condition
      Error in `read_multi()`:
      ! `paths` must be a named character vector.

# compute_charging_diffs errors when condition column missing

    Code
      compute_charging_diffs(charging, numerator = "a", denominator = "b")
    Condition
      Error in `compute_charging_diffs()`:
      ! Column "condition" not found in `charging_data`.

