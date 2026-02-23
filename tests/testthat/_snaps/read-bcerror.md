# compute_bcerror_delta errors on invalid delta expression

    Code
      compute_bcerror_delta(df, delta = wt + mut)
    Condition
      Error in `compute_bcerror_delta()`:
      ! `delta` must be an expression of the form `lhs - rhs`.

# compute_bcerror_delta errors on missing condition level

    Code
      compute_bcerror_delta(df, delta = wt - missing)
    Condition
      Error in `compute_bcerror_delta()`:
      ! Level "missing" not found in column condition.

