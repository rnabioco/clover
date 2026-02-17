# Propagate standard error for a difference.

Compute the standard error of `a - b` assuming independent errors.

## Usage

``` r
propagate_error_diff(se_a, se_b)
```

## Arguments

- se_a:

  Numeric vector of standard errors for the first value.

- se_b:

  Numeric vector of standard errors for the second value.

## Value

Numeric vector of propagated standard errors.

## Examples

``` r
propagate_error_diff(1, 0.5)
#> [1] 1.118034
```
