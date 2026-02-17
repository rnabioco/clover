# Propagate standard error for a ratio.

Compute the standard error of `a / b` using first-order error
propagation.

## Usage

``` r
propagate_error_ratio(a, b, se_a, se_b)
```

## Arguments

- a:

  Numeric vector of numerator values.

- b:

  Numeric vector of denominator values.

- se_a:

  Numeric vector of standard errors for `a`.

- se_b:

  Numeric vector of standard errors for `b`.

## Value

Numeric vector of propagated standard errors.

## Examples

``` r
propagate_error_ratio(10, 5, 1, 0.5)
#> [1] 0.2828427
```
