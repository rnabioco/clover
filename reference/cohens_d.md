# Calculate Cohen's d effect size.

Compute the standardized mean difference between two groups using the
pooled standard deviation.

## Usage

``` r
cohens_d(x, y)
```

## Arguments

- x:

  Numeric vector for the first group.

- y:

  Numeric vector for the second group.

## Value

A single numeric value for the effect size.

## Examples

``` r
cohens_d(rnorm(20, mean = 5), rnorm(20, mean = 3))
#> [1] 1.58177
```
