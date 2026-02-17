# Calculate log2 fold change with pseudocount.

Calculate log2 fold change with pseudocount.

## Usage

``` r
calc_fold_change(x, y, pseudocount = 1)
```

## Arguments

- x:

  Numeric vector of numerator values.

- y:

  Numeric vector of denominator values.

- pseudocount:

  Pseudocount added to both numerator and denominator before computing
  the ratio. Default `1`.

## Value

Numeric vector of log2 fold changes.

## Examples

``` r
calc_fold_change(10, 5)
#> [1] 0.8744691
calc_fold_change(c(0, 10), c(10, 0), pseudocount = 1)
#> [1] -3.459432  3.459432
```
